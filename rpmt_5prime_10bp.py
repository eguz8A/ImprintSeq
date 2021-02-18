
import sys

read_id = ""
read_seq = ""
read_3 = ""
read_qual = ""

out_file = open(sys.argv[1].replace(".fq", "_rpmt.fq"), 'w')

with open(sys.argv[1]) as infile:
    line_counter = 0
    for line in infile:
        if line_counter%4 == 0 and line.startswith("@"):
            read_id = line.strip()
        elif line_counter%4 == 1:
            read_seq = line.strip()
        elif line_counter%4 == 2 and line.startswith("+"):
            read_3 = line.strip()
        else:
            read_qual = line.strip()

            out_file.write(read_id + ":" + read_seq[:10] + "\n")
            out_file.write(read_seq[10:] + "\n")
            out_file.write(read_3 + "\n")
            out_file.write(read_qual[10:] + "\n")

        line_counter += 1

out_file.close()