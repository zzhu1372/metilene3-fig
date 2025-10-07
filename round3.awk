#!/usr/bin/awk -f
BEGIN {
    FS = OFS = "\t"
}

# Print the first line (header) unchanged
NR == 1 {
    print
    next
}

{
    # Print the first 3 columns as-is
    printf "%s\t%s\t%s", $1, $2, $3

    # Format remaining columns to 3 decimal digits
    for (i = 4; i <= NF; i++) {
        if ($i == ".")
            printf "\t."
        else
            printf "\t%.3f", $i
    }

    # End line
    printf "\n"
}
