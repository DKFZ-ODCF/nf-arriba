#!/usr/bin/awk -f

BEGIN { FS="\t" } # Sophia files are tab-separated

NR==1 { for (i=1; i<=NF; ++i) col[$i]=i } # get column names

$col["source1"] ~ /\|/ && $col["source2"] ~ /\|/ { # skip lines without direction info (i.e., lacking the pipe symbol)

	# extract coordinates
	position1 = $col["#chrom1"] ":" $col["end1"]
	position2 = $col["chrom2"] ":" $col["end2"]

	# convert SOPHIA encoding of breakpoint orientations to Arriba:
	# pipe at the beginning means downstream; anywhere else means upstream
	direction1 = ($col["source1"] ~ /^\|/) ? "downstream" : "upstream"
	direction2 = ($col["source2"] ~ /^\|/) ? "downstream" : "upstream"

	print position1 "\t" position2 "\t" direction1 "\t" direction2
}

