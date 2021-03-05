import sys, os

if __name__ == "__main__":

	with sys.stdin as fin, sys.stdout as fout:

		addAAtext = False
		for line in fin:
			if line.startswith("#"):
				if line.startswith("##INFO"):
					if not addAAtext:
						fout.write("##INFO=<ID=AA,Number=A,Type=String,Description=\"Ancestral state based on chimpanzee seq\">"+"\n")
						addAAtext = True
				fout.write(line)
				
			else:
				tmp_line = line.strip().split()
				tmp_line[7] = "AA=0;" + tmp_line[7]
				fout.write("\t".join(tmp_line) + "\n")
