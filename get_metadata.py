# Author: Leandro de Mattos Pereira
# 06 Nov 2024
# Team: CNP Team - Leam Lab
#
# Updated: 15 Sep 2025
# - Resumes automatically by reading what already exists in the output file (and optionally from --resume).
# - At the END of the run, the filtered file is REBUILT FROM SCRATCH and OVERWRITES the current one:
#   (filtered_<output>.tsv), with the criterion: BiomeDistribution != 'Unknown' AND Latitude/Longitude != 'Unknown'.
# - NEW: --retry-only (heal/complete only the existing output) and --retry-missing (heal after the normal flow).
# - NEW: --retry-columns to limit which fields to heal; --retry-changed-only to log only improved rows.
#
# Notes:
# - <input_assemblies.tsv> must have Assembly Accessions in the 1st column; by default we assume a HEADER.
#   If your IDs file does NOT have a header, use --no-input-header.
# - Entrez API key (optional): export NCBI_API_KEY="YOUR_KEY"
# - Email can come from --email or env ENTREZ_EMAIL (fallback to default in code).
# - Requirements: NCBI CLI 'datasets' and 'dataformat' in PATH; ete3; biopython.

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional, Set, Tuple, Dict, List

from ete3 import NCBITaxa
from Bio import Entrez
import xml.etree.ElementTree as ElementTree

# ======================
# Entrez configuration
# ======================

def configure_entrez(email: Optional[str] = None):
  default_email = 'mattoslmp@gmail.com'
  Entrez.email = email or os.getenv('ENTREZ_EMAIL') or default_email
  api_key = os.getenv('NCBI_API_KEY')
  if api_key:
    Entrez.api_key = api_key

# ======================
# File utilities
# ======================

def file_has_content(p: Path) -> bool:
  try:
    return p.exists() and p.stat().st_size > 0
  except Exception:
    return False

def read_processed_accessions(metadata_path: Path) -> Set[str]:
  """
  Reads a previously generated metadata TSV (same layout) and returns
  the set of 'Assembly Accession' already present (1st column).
  Ignores header and blank lines.
  """
  processed = set()
  if not file_has_content(metadata_path):
    return processed
  try:
    with metadata_path.open('r', encoding='utf-8', errors='replace') as f:
      first = True
      for line in f:
        if not line.strip():
          continue
        if first:
          first = False
          if line.lower().startswith("assembly accession"):
            continue
        parts = line.rstrip('\n').split('\t')
        if parts:
          processed.add(parts[0])
  except Exception as e:
    print(f"[Warn] Could not read '{metadata_path}': {e}", file=sys.stderr)
  return processed

def is_header_line(line: str) -> bool:
  l = line.strip().lower()
  return (
    l.startswith("assembly accession")
    or ("assembly" in l and "accession" in l)
    or ("organism" in l and "tax id" in l)
  )

# ======================
# Taxonomy (ETE3)
# ======================

ncbi = NCBITaxa()

def get_lineage(tax_id: int) -> str:
  try:
    lineage = ncbi.get_lineage(tax_id)
    names = ncbi.get_taxid_translator(lineage)
    lineage_str = "; ".join([names[t] for t in lineage])
    return lineage_str
  except Exception as e:
    print(f"[Error] Failed to get lineage for Tax ID {tax_id}: {e}")
    return "Unknown"

# ======================
# Biome / Coordinates
# ======================

def parse_latitude_longitude(lat_lon_str: Optional[str]) -> Tuple[Optional[float], Optional[float]]:
  if lat_lon_str and isinstance(lat_lon_str, str):
    try:
      # Typical format: "12.34 N 56.78 W"
      parts = lat_lon_str.strip().replace(',', '.').split()
      if len(parts) >= 4:
        latitude_value = float(parts[0])
        latitude_direction = parts[1].upper()
        longitude_value = float(parts[2])
        longitude_direction = parts[3].upper()

        if latitude_direction == 'S':
          latitude_value = -latitude_value
        if longitude_direction == 'W':
          longitude_value = -longitude_value

        return latitude_value, longitude_value
    except Exception as e:
      print(f"[Warn] Could not parse lat_lon '{lat_lon_str}': {e}")
  return None, None

def categorize_biome(biome_description: Optional[str]) -> str:
  if biome_description is None:
    return 'Unknown'
  biome_description = biome_description.lower()
  # Coarse biome buckets based on common GOLD-like descriptors
  if any(term in biome_description for term in ['soil', 'forest', 'desert', 'savanna']):
    return 'Terrestrial'
  elif any(term in biome_description for term in ['marine', 'sea', 'ocean', 'coastal']):
    return 'Marine'
  elif any(term in biome_description for term in ['lake', 'freshwater', 'river', 'pond']):
    return 'Freshwater'
  elif any(term in biome_description for term in ['waste', 'wastewater', 'sewage']):
    return 'Wastewater'
  elif any(term in biome_description for term in ['host', 'symbiont', 'root', 'nodule']):
    return 'Host-Associated'
  elif 'hypersaline' in biome_description:
    return 'Extreme - Hypersaline'
  elif 'hot spring' in biome_description or 'thermal' in biome_description:
    return 'Extreme - Thermal'
  elif 'acidic' in biome_description or 'alkaline' in biome_description:
    return 'Extreme - Acidic/Alkaline'
  elif any(term in biome_description for term in ['reef', 'coral']):
    return 'Reef'
  elif 'environmental sample' in biome_description:
    return 'Environmental Sample'
  else:
    return 'Other'

# ======================
# Entrez / BioSample
# ======================

# Counters for this run
biome_count = 0
lat_lon_count = 0

def fetch_biosample_metadata(assembly_accession: str, delay: float = 0.35):
  global biome_count, lat_lon_count
  try:
    search_handle = Entrez.esearch(db='assembly', term=f'{assembly_accession}[Assembly Accession]')
    search_record = Entrez.read(search_handle)
    search_handle.close()
    time.sleep(delay)

    if search_record.get('IdList'):
      assembly_uid = search_record['IdList'][0]
      link_handle = Entrez.elink(dbfrom='assembly', id=assembly_uid, db='biosample')
      link_records = Entrez.read(link_handle)
      link_handle.close()
      time.sleep(delay)

      if link_records and link_records[0].get('LinkSetDb'):
        biosample_uid = link_records[0]['LinkSetDb'][0]['Link'][0]['Id']
        fetch_handle = Entrez.efetch(db='biosample', id=biosample_uid, rettype='xml')
        xml_data = fetch_handle.read()
        fetch_handle.close()
        time.sleep(delay)

        root = ElementTree.fromstring(xml_data)
        latitude_longitude = None
        environmental_sample_flag = None
        isolation_source = None
        geographic_location_name = None

        # Extract attributes from XML
        attributes = root.findall('.//Attribute')
        for attribute in attributes:
          attribute_name = attribute.get('attribute_name', '').lower()
          if attribute_name == 'lat_lon':
            latitude_longitude = attribute.text
          elif attribute_name == 'environmental_sample':
            environmental_sample_flag = attribute.text
          elif attribute_name == 'isolation_source':
            isolation_source = attribute.text
          elif attribute_name == 'geo_loc_name':
            geographic_location_name = attribute.text

        if environmental_sample_flag and str(environmental_sample_flag).lower() == 'true':
          environmental_sample_flag = 'Environmental sample'

        biome_description = isolation_source or environmental_sample_flag
        latitude, longitude = parse_latitude_longitude(latitude_longitude)

        if biome_description:
          biome_count += 1
        if latitude is not None and longitude is not None:
          lat_lon_count += 1

        return {
          'Location': geographic_location_name or "Unknown",
          'BiomeDistribution': categorize_biome(biome_description),
          'Latitude': latitude if latitude is not None else "Unknown",
          'Longitude': longitude if longitude is not None else "Unknown"
        }
  except Exception as e:
    print(f"[Error] Failed to fetch BioSample for {assembly_accession}: {e}")
  return {
    'Location': "Unknown",
    'BiomeDistribution': "Unknown",
    'Latitude': "Unknown",
    'Longitude': "Unknown"
  }

# ======================
# NCBI Datasets / dataformat
# ======================

def run_ncbi_datasets_summary(assembly_id: str) -> Tuple[int, str, str]:
  """
  Runs datasets + dataformat and returns (returncode, stdout, stderr)
  """
  command = (
    f"datasets summary genome accession {assembly_id} --as-json-lines | "
    f"dataformat tsv genome --fields "
    f"accession,organism-name,organism-common-name,organism-tax-id,"
    f"assminfo-level,assminfo-bioproject,assminfo-biosample-accession,"
    f"assmstats-gc-percent,assmstats-total-sequence-len,assminfo-sequencing-tech,"
    f"assminfo-release-date,assminfo-biosample-collection-date,assminfo-biosample-description-title"
  )
  result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  return result.returncode, result.stdout.decode('utf-8', errors='replace'), result.stderr.decode('utf-8', errors='replace')

# ======================
# Header / writing
# ======================

HEADER = (
  "Assembly Accession\tOrganism Name\tOrganism Common Name\tOrganism Tax ID\tLineage\t"
  "Assembly Level\tBioProject Accession\tBioSample Accession\tGC Percent\t"
  "Total Sequence Length\tSequencing Technology\tRelease Date\tCollection Date\t"
  "BioSample Description\tLocation\tBiomeDistribution\tLatitude\tLongitude\n"
)

def ensure_header(path: Path):
  """Writes header only if the file is empty/non-existent."""
  if not file_has_content(path):
    with path.open('w', encoding='utf-8') as f:
      f.write(HEADER)

def append_line(path: Path, line: str):
  with path.open('a', encoding='utf-8') as f:
    f.write(line + "\n")

def rebuild_filtered(full_output_path: Path, filtered_output_path: Path):
  """
  Rebuilds the filtered file FROM SCRATCH and OVERWRITES the existing one.
  Criterion: BiomeDistribution != 'Unknown' AND Latitude/Longitude != 'Unknown'.
  """
  tmp_path = filtered_output_path.with_suffix(filtered_output_path.suffix + ".tmp")
  with full_output_path.open('r', encoding='utf-8', errors='replace') as fin, \
       tmp_path.open('w', encoding='utf-8') as fout:
    # header
    fout.write(HEADER)
    first = True
    for line in fin:
      if not line.strip():
        continue
      if first:
        first = False
        if line.lower().startswith("assembly accession"):
          continue
      parts = line.rstrip('\n').split('\t')
      if len(parts) < 18:
        continue
      biome = parts[-3].strip()
      lat = parts[-2].strip()
      lon = parts[-1].strip()
      if biome != "Unknown" and lat != "Unknown" and lon != "Unknown":
        fout.write(line if line.endswith('\n') else line + '\n')
  os.replace(tmp_path, filtered_output_path)

# ======================
# Retry missing fields
# ======================

# Column indices according to HEADER:
IDX_ACCESSION = 0
IDX_TAXID = 3
IDX_LINEAGE = 4
IDX_LOCATION = -4
IDX_BIOME = -3
IDX_LAT = -2
IDX_LON = -1

def parse_retry_columns(s: Optional[str]) -> Set[str]:
  """
  Converts the CSV string from --retry-columns to a normalized set.
  Accepted values: biome, latlon, location, lineage
  """
  if not s:
    return {"biome", "latlon", "location", "lineage"}
  allowed = {"biome", "latlon", "location", "lineage"}
  out = set()
  for tok in s.split(','):
    t = tok.strip().lower()
    if t in allowed:
      out.add(t)
  return out or {"biome", "latlon", "location", "lineage"}

def needs_retry(parts: List[str], cols: Set[str]) -> bool:
  try:
    lineage_bad = (parts[IDX_LINEAGE].strip() == "Unknown")
    loc_bad = (parts[IDX_LOCATION].strip() == "" or parts[IDX_LOCATION].strip() == "Unknown")
    biome_bad = (parts[IDX_BIOME].strip() == "" or parts[IDX_BIOME].strip() == "Unknown")
    lat_bad = (parts[IDX_LAT].strip() == "" or parts[IDX_LAT].strip() == "Unknown")
    lon_bad = (parts[IDX_LON].strip() == "" or parts[IDX_LON].strip() == "Unknown")
    latlon_bad = lat_bad or lon_bad
    return (
      ("lineage" in cols and lineage_bad)
      or ("location" in cols and loc_bad)
      or ("biome" in cols and biome_bad)
      or ("latlon" in cols and latlon_bad)
    )
  except Exception:
    return False

def merge_retry(parts: List[str], tax_id_str: str, meta: Dict[str, str], cols: Set[str]) -> Tuple[List[str], List[str]]:
  """
  Updates fields allowed by 'cols'. Returns (updated_parts, improved_fields_list)
  """
  improved: List[str] = []

  # Lineage
  if "lineage" in cols and parts[IDX_LINEAGE].strip() == "Unknown" and tax_id_str.isdigit():
    lineage_new = get_lineage(int(tax_id_str))
    if lineage_new and lineage_new != "Unknown":
      parts[IDX_LINEAGE] = lineage_new
      improved.append("lineage")

  # BioSample-derived fields
  if meta:
    if "location" in cols and (parts[IDX_LOCATION].strip() in ("", "Unknown")) and meta.get("Location"):
      parts[IDX_LOCATION] = str(meta["Location"])
      improved.append("location")
    if "biome" in cols and (parts[IDX_BIOME].strip() in ("", "Unknown")) and meta.get("BiomeDistribution"):
      parts[IDX_BIOME] = str(meta["BiomeDistribution"])
      improved.append("biome")
    if "latlon" in cols:
      if parts[IDX_LAT].strip() in ("", "Unknown") and str(meta.get("Latitude", "")) not in ("", "Unknown"):
        parts[IDX_LAT] = str(meta["Latitude"])
        if "latlon" not in improved:
          improved.append("latlon")
      if parts[IDX_LON].strip() in ("", "Unknown") and str(meta.get("Longitude", "")) not in ("", "Unknown"):
        parts[IDX_LON] = str(meta["Longitude"])
        if "latlon" not in improved:
          improved.append("latlon")

  return parts, improved

def retry_fill_missing(output_path: Path, delay: float, cols: Set[str], log_changed_only: bool):
  """
  Rewrites 'output_path' trying to fill missing fields per row.
  Uses Entrez/ETE3 only (does NOT re-run datasets/dataformat).
  """
  tmp_path = output_path.with_suffix(output_path.suffix + ".retry.tmp")
  updated = 0
  total_checked = 0

  with output_path.open('r', encoding='utf-8', errors='replace') as fin, \
       tmp_path.open('w', encoding='utf-8') as fout:

    first = True
    for raw in fin:
      line = raw.rstrip('\n')
      if not line.strip():
        continue

      if first:
        first = False
        if line.lower().startswith("assembly accession"):
          fout.write(line + "\n")
          continue

      parts = line.split('\t')
      if len(parts) < 18:
        fout.write(line + "\n")
        continue

      total_checked += 1

      if not needs_retry(parts, cols):
        fout.write(line + "\n")
        continue

      accession = parts[IDX_ACCESSION].strip()
      tax_id_str = parts[IDX_TAXID].strip()

      meta = fetch_biosample_metadata(accession, delay=delay) if ({"biome","latlon","location"} & cols) else {}
      parts_new, improved = merge_retry(parts, tax_id_str, meta, cols)

      if improved:
        updated += 1
        if log_changed_only:
          print(f"[retry] {accession} improved: {','.join(improved)}")
      else:
        if not log_changed_only:
          print(f"[retry] {accession} no improvement")

      fout.write("\t".join(parts_new) + "\n")

  os.replace(tmp_path, output_path)
  print(f"[retry] Checked rows: {total_checked}")
  print(f"[retry] Updated rows: {updated}")

# ======================
# Main pipeline
# ======================

def main():
  parser = argparse.ArgumentParser(
    description="Fetches NCBI assembly metadata, appends to a single TSV, and rebuilds a filtered TSV at the end."
  )
  parser.add_argument("input_file", help="Input TSV with Assembly Accessions in the 1st column.")
  parser.add_argument("output_file", help="Sole metadata TSV (will be appended).")
  parser.add_argument("--resume", help="Additional pre-existing TSV used to skip already processed accessions (optional).", default=None)
  parser.add_argument("--no-input-header", help="Use if the input IDs file has NO header.", action="store_true")
  parser.add_argument("--email", help="Entrez email (takes precedence over env ENTREZ_EMAIL).", default=None)
  parser.add_argument("--delay", help="Delay (seconds) between Entrez calls.", type=float, default=0.35)

  # Healing flags
  parser.add_argument("--retry-missing", help="After processing, try to fill missing fields in output_file (Biome/Lat/Lon/Location/Lineage).", action="store_true")
  parser.add_argument("--retry-only", help="Do not process new IDs; only heal missing fields in output_file and exit.", action="store_true")
  parser.add_argument("--retry-columns", help="Limit which fields to heal (CSV): biome,latlon,location,lineage. Default: all.", default=None)
  parser.add_argument("--retry-changed-only", help="During retry, log only rows that were improved.", action="store_true")

  args = parser.parse_args()

  configure_entrez(email=args.email)

  input_path = Path(args.input_file)
  output_path = Path(args.output_file)
  filtered_output_path = output_path.parent / f"filtered_{output_path.name}"

  # Ensure header on main output
  ensure_header(output_path)

  # Columns to heal
  retry_cols = parse_retry_columns(args.retry_columns)

  # Retry-only mode
  if args.retry_only:
    retry_fill_missing(output_path, delay=args.delay, cols=retry_cols, log_changed_only=args.retry_changed_only)
    rebuild_filtered(output_path, filtered_output_path)
    print("[retry-only] Done. Filtered file rebuilt.")
    return

  # Processed accession set (output + optional resume)
  processed: Set[str] = set()
  processed |= read_processed_accessions(output_path)
  if args.resume:
    processed |= read_processed_accessions(Path(args.resume))

  # Stats
  skipped_already = 0
  processed_now = 0

  if not input_path.exists():
    print(f"[Error] Input file '{input_path}' not found.", file=sys.stderr)
    sys.exit(1)

  with input_path.open('r', encoding='utf-8', errors='replace') as in_file:
    first_line = True
    for raw in in_file:
      line = raw.rstrip('\n')
      if not line.strip():
        continue

      if first_line:
        first_line = False
        if not args.no_input_header and is_header_line(line):
          continue  # skip header

      # 1st column = accession (TSV or whitespace)
      assembly_id = (line.split('\t')[0] if '\t' in line else line.split()[0]).strip()
      if not assembly_id:
        continue

      if assembly_id in processed:
        skipped_already += 1
        continue

      # datasets/dataformat
      ret, stdout, stderr = run_ncbi_datasets_summary(assembly_id)
      if ret != 0:
        print(f"[Error] Processing {assembly_id}: {stderr}")
        continue

      output_lines = stdout.splitlines()
      if output_lines and output_lines[0].lower().startswith("assembly accession"):
        output_lines = output_lines[1:]
      if not output_lines:
        print(f"[Warn] No lines returned for {assembly_id}.")
        continue

      wrote_any_for_assembly = False

      for output_line in output_lines:
        if not output_line.strip():
          continue

        fields = output_line.split('\t')
        if len(fields) < 13:
          print(f"[Warn] Unexpected line for {assembly_id}: {output_line}")
          continue

        tax_id = fields[3] if len(fields) > 3 else ""
        lineage = get_lineage(int(tax_id)) if tax_id.isdigit() else "Unknown"

        # Additional metadata (BioSample)
        metadata = fetch_biosample_metadata(assembly_id, delay=args.delay)

        # Enriched line
        enriched_line = (
          "\t".join(fields[:4]) + f"\t{lineage}\t" +
          "\t".join(fields[4:]) + f"\t{metadata['Location']}\t{metadata['BiomeDistribution']}\t"
          f"{metadata['Latitude']}\t{metadata['Longitude']}"
        )

        append_line(output_path, enriched_line)
        wrote_any_for_assembly = True

      if wrote_any_for_assembly:
        processed.add(assembly_id)
        processed_now += 1

  # Optional healing after normal processing
  if args.retry_missing:
    retry_fill_missing(output_path, delay=args.delay, cols=retry_cols, log_changed_only=args.retry_changed_only)

  # Rebuild filtered file
  rebuild_filtered(output_path, filtered_output_path)

  # Summary
  print(f"Already present & skipped: {skipped_already}")
  print(f"Processed in this run: {processed_now}")
  print(f"With BiomeDistribution filled in this run: {biome_count}")
  print(f"With Latitude/Longitude filled in this run: {lat_lon_count}")
  print(f"Main output (appended): {output_path}")
  print(f"Filtered rebuilt (overwritten): {filtered_output_path}")

if __name__ == "__main__":
  main()
