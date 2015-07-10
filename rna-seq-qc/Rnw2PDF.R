require(knitr)

args = commandArgs(TRUE)

cwd = args[1]

## set variable for Rnw code chunk
opts_knit$set(report.dir = cwd)  #report dir
# opts_knit$get()

# infile_Rnw = "/home/kilpert/git/rna-seq-qc/rna-seq-qc/Report.Rnw"
# outfile_tex = "/data/processing/kilpert/test/rna-seq-qc/Ausma/PE_mm10_FULL_hisat_trim/project_report/Report.tex"
# outfile_pdf = "/data/processing/kilpert/test/rna-seq-qc/Ausma/PE_mm10_FULL_hisat_trim/project_report/Report.pdf"

infile_Rnw = args[2]
outfile_tex = file.path(cwd, "Report.tex")
outfile_pdf = file.path(cwd, "Report.pdf")

knit(infile_Rnw, output=outfile_tex)    ## build .tex file
system( sprintf("pdflatex -output-directory=%s %s ", cwd, outfile_tex) )  ##  make .pdf
