import xml.etree.ElementTree as ET

def extract_entry_by_name(xml_file, target_name, output_file):
    ns = {"ns": "https://uniprot.org/uniprot"}
    context = ET.iterparse(xml_file, events=("start", "end"))
    context = iter(context)
    _, root = next(context)

    found = False
    for event, elem in context:
        if event == "end" and elem.tag == "{https://uniprot.org/uniprot}entry":
            entry_name = elem.findtext("ns:name", namespaces=ns)
            if entry_name == target_name:
                print(f"✅ Found entry: {entry_name}")
                found = True

                # Wrap it in a <uniprot> root for a valid XML document
                wrapper = ET.Element("uniprot", xmlns=ns["ns"])
                wrapper.append(elem)

                tree = ET.ElementTree(wrapper)
                tree.write(output_file, encoding="utf-8", xml_declaration=True)
                break
            else:
                root.clear()

    if not found:
        print(f"❌ Entry {target_name} not found.")

#extract_entry_by_name("../0-download-inputs/data-files/uniprot_sprot.xml", "ACEA_YEAST", "acea_yeast_entry.xml")
extract_entry_by_name("../0-download-inputs/data-files/uniprot_sprot.xml", "CHO2_YEAST", "cho2_yeast_entry.xml")
