

{ “cells”: \[ { “cell_type”: “markdown”, “metadata”: {}, “source”: \[
“—”, “title: "Class 16"”, “author: "Jiayi Zhou (PID:A17856751)"”, “date:
"2025-11-20"”, “output: pdf_document”, “—”, “”, “\## Read the BLAST TSV
and set column names”, “”,
“`{r}\n",         "b <- read.table(\"mm-second.x.zebrafish.tsv\",\n",         "                sep = \"\\t\",\n",         "                header = FALSE,\n",         "                stringsAsFactors = FALSE)\n",         "`”,
“”,
“`{r}\n",         "colnames(b) <- c(\n",         "  \"qseqid\", \"sseqid\", \"pident\", \"length\",\n",         "  \"mismatch\", \"gapopen\", \"qstart\", \"qend\",\n",         "  \"sstart\", \"send\", \"evalue\", \"bitscore\"\n",         ")\n",         "`”,
“”, “`{r}\n",         "head(b)\n",         "str(b)\n",         "`”, “\##
Make a histogram of bitscore”, “”,
“`{r}\n",         "hist(b$bitscore,\n",         "     breaks = 30,\n",         "     main = \"Distribution of BLAST bitscores\",\n",         "     xlab = \"Bitscore\")\n",         "`”,
“”, “Larger bitscores are better.”, “”,
“`{r}\n",         "adj_identity <- b$pident * (b$qend - b$qstart)\n",         "\n",         "plot(adj_identity, b$bitscore,\n",         "     xlab = \"pident × alignment length (qend - qstart)\",\n",         "     ylab = \"bitscore\",\n",         "     main = \"Adjusted identity vs bitscore\")\n",         "`”,
“”, “\## Using ggplot2”, “”,
“`{r}\n",         "library(ggplot2)\n",         "`”, “”,
“`{r}\n",         "ggplot(b, aes(x = pident, y = bitscore)) +\n",         "  geom_point(alpha = 0.1) +\n",         "  labs(title = \"Percent identity vs bitscore\")\n",         "`”,
“”,
“`{r}\n",         "ggplot(b, aes(x = pident * (qend - qstart), y = bitscore)) +\n",         "  geom_point(alpha = 0.1) +\n",         "  geom_smooth() +\n",         "  labs(title = \"pident × alignment length vs bitscore\",\n",         "       x = \"pident × (qend - qstart)\",\n",         "       y = \"bitscore\")\n",         "`”,
“Bitscore increases with both identity and alignment length.”, “Two
alignments with the same percent identity can have different bitscores
if one covers more residues.”, “So bitscore is only somewhat related to
percent identity alone; it’s more tightly linked to “how many residues
are matching over how long a region.”“,”“,”\## Knit the R
Markdown“,”“,”Run
`scp -i ~/Downloads/keyjz.pem -r ubuntu@ec2-44-252-77-168.us-west-2.compute.amazonaws.com:~/work .`
in terminal.“,”“,”\> `-r` purpose“,”“,”- Means recursive.“,”“,”- It
copies directories and all their contents (subdirectories, files,
etc.).“,”“,”- Without `-r`, scp will fail on directories.“,”“,”\> `*`
Pupose“,” “,”- `~/work/*` means “all files and directories inside
~/work”.“,”“,”- Copies everything in `~/work` from the remote machine
into current local directory.“,”“,”\## Using rsync“,”“,”Run
`scp -i ~/Downloads/keyjz.pem -r ubuntu@ec2-44-252-77-168.us-west-2.compute.amazonaws.com:~/work .`
in terminal.“,”“,”- **`-a`**: Archive mode – copy recursively and
preserve permissions, timestamps, and other metadata so the directory is
cloned as faithfully as possible.“,”“,”- **`-z`**: Compress file data
during transfer to reduce network bandwidth usage and often speed up
transfers.“,”“,”- **`-P`**: Show progress for each file and keep
partially transferred files so interrupted transfers can resume.“,”“,”-
**`--exclude`**: Skip copying any files or directories whose names match
the given pattern (e.g. `--exclude=\"*.psq\"`).” \] } \], “metadata”: {
“kernelspec”: { “display_name”: “Python 3”, “language”: “python”,
“name”: “python3” } }, “nbformat”: 4, “nbformat_minor”: 4 }
