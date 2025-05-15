import xml.etree.ElementTree as ET
import pandas as pd
import argparse
import re
from goatools.obo_parser import GODag

# sloppy hardcoding but I just want this to work
go_dag = GODag("0-download-inputs/data-files/go-basic.obo")

def parse_entry(entry, target_organism):

    ns = {"ns": "https://uniprot.org/uniprot"}
    organism_name = entry.findtext("ns:organism/ns:name[@type='scientific']", namespaces=ns)

    if organism_name != target_organism:
        return None

    accession = entry.findtext("ns:accession", namespaces=ns)
    entry_name = entry.findtext("ns:name", namespaces=ns)
    protein_name = entry.findtext("ns:protein/ns:recommendedName/ns:fullName", namespaces=ns)
    sequence = entry.findtext("ns:sequence", namespaces=ns)

    # extract GO terms
    go_terms = [
        db_ref.attrib["id"]
        for db_ref in entry.findall("ns:dbReference", namespaces=ns)
        if db_ref.attrib.get("type") == "GO"
    ]

    # and map them to human-readable names with GOATOOLS
    go_mapped = [
        f"{go_id}: {go_dag[go_id].name}" if go_id in go_dag else go_id
        for go_id in go_terms
    ]

    # extract FT-based PTMs
    ptm_keywords = {"modified residue", "lipid moiety-binding region", "glycosylation site"}
    ft_ptms = []
    for ft in entry.findall("ns:feature", namespaces=ns):
        ft_type = ft.attrib.get("type")
        if ft_type in ptm_keywords:
            desc = ft.attrib.get("description", "")
            # get position (either <position> or range <begin>/<end>)
            loc_elem = ft.find("ns:location/ns:position", namespaces=ns)
            if loc_elem is not None:
                position = loc_elem.attrib.get("position")
            else:
                # handle range-based features as fallback
                begin = ft.find("ns:location/ns:begin", namespaces=ns)
                end = ft.find("ns:location/ns:end", namespaces=ns)
                if begin is not None and end is not None:
                    position = f"{begin.attrib.get('position')}-{end.attrib.get('position')}"
                else:
                    position = "?"

            # format output
            ft_ptms.append(f"{desc} at {position}")

    # extract function comments
    function_texts = []
    for comment in entry.findall("ns:comment[@type='function']", namespaces=ns):
        for text_elem in comment.findall("ns:text", namespaces=ns):
            text = text_elem.text
            if text:
                # remove citation notes like "{ECO:...}"
                cleaned = re.sub(r"\{.*?\}", "", text).strip().rstrip(".")
                function_texts.append(cleaned)

    # extract subcellular location comments
    subcell_locs = []
    for comment in entry.findall("ns:comment[@type='subcellular location']", namespaces=ns):
        for subloc in comment.findall("ns:subcellularLocation", namespaces=ns):
            loc = subloc.findtext("ns:location", namespaces=ns)
            topology = subloc.findtext("ns:topology", namespaces=ns)
            orientation = subloc.findtext("ns:orientation", namespaces=ns)
            parts = [loc]
            if topology:
                parts.append(f"topology: {topology}")
            if orientation:
                parts.append(f"orientation: {orientation}")
            joined = " (" + "; ".join(parts[1:]) + ")" if len(parts) > 1 else ""
            if loc:
                subcell_locs.append(loc + joined)

    return {
        "EntryName": entry_name,
        "PrimaryAccession": accession,
        "ProteinName": protein_name,
        "Organism": organism_name,
        "Sequence": sequence,
        "GO_terms": ";".join(go_terms),
        "GO_terms_human_readable": ";".join(go_mapped),
        "parsed_PTMs": ";".join(ft_ptms),
        "parsed_functions": ";".join(function_texts),
        "localization_keywords": ";".join(subcell_locs),
    }


def parse_uniprot_xml(xml_file, organism_filter):
    ns = {"ns": "https://uniprot.org/uniprot"}
    entries = []

    context = ET.iterparse(xml_file, events=("start", "end"))
    context = iter(context)
    _, root = next(context)

    for event, elem in context:
        if event == "end" and elem.tag == "{https://uniprot.org/uniprot}entry":
            entry_name = elem.findtext("ns:name", namespaces=ns)
            organism = elem.findtext("ns:organism/ns:name[@type='scientific']", namespaces=ns)

            #if entry_name == "ACEA_YEAST":
            #    print(f"✅ Found target entry: {entry_name}")
            #    print(f"Organism: {organism}")
            #    for ft in elem.findall("ns:feature", namespaces=ns):
            #        if ft.attrib.get("type") == "modified residue":
            #            position = ft.find("ns:location/ns:position", namespaces=ns)
            #            desc = ft.findtext("ns:description", default="", namespaces=ns)
            #            pos_val = position.attrib["position"] if position is not None else "?"
            #            print(f"Feature: modified residue at {pos_val}")
            #            print(f"Description: '{desc}'")

            record = parse_entry(elem, organism_filter)
            if record:
                entries.append(record)

            root.clear()

    return pd.DataFrame(entries)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True)
    parser.add_argument("--output_file", type=str, required=True)
    parser.add_argument("--organism", type=str, required=True)
    args = parser.parse_args()

    df = parse_uniprot_xml(args.input_file, args.organism)
    df.to_csv(args.output_file, index=False)


if __name__ == "__main__":
    main()
