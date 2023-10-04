# Script used to modify gene names in synthetic Chr11 GFF3 file to distinguish 
# from wild-type Chr11 genes.

# Process the file using awk
awk -F '\t' 'BEGIN { OFS = "\t" }

# skip first 3 lines
NR <= 3 { print; next }

{
    # split attributes field into individual attributes
    split($9, attributes, ";")
    new_attributes = ""
    for (i = 1; i <= length(attributes); i++) {
    attr = attributes[i]

    # check if the attribute starts with "ID=" or "Name="
    if (match(attr, /^ID=/) || match(attr, /^Name=/)) {

        # extract the attribute name and value
        attr_name = substr(attr, 1, RSTART + RLENGTH - 2)
        attr_value = substr(attr, RLENGTH + 1)

        # check if the attribute name is "ID" or "Name" and if it has a value
        if ((attr_name == "ID" || attr_name == "Name") && length(attr_value) > 0) {

        # check if "gene:" or "transcript:" prefix is already present
        if (!match(attr_value, /^(gene|transcript):/)) {
            attr_value = "x." attr_value
        } else {
            attr_value = substr(attr_value, 1, RLENGTH) "x." substr(attr_value, RLENGTH + 1)
        }
        }

        # reconstruct the modified attribute
        modified_attr = attr_name "=" attr_value
    } else {
        modified_attr = attr
    }
    new_attributes = new_attributes ";" modified_attr
    }

    # replace the attributes field with the modified attributes
    $9 = substr(new_attributes, 2)
    print
}
' tmp-syn.gff > modified_syn_chr11.gff
