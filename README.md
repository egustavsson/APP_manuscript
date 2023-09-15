
<!-- README.md is generated from README.Rmd. Please edit that file -->
# APP transcriptomic manuscript

## Background

This repository contains all the scripts used to in the APP manuscript post transcript calling.

### Processing of long-read data

We have previosly developed two `snakemake` pipelines for automated Analysis of Iso-Seq data which we have used here:

[APTARS (Analysis of PacBio TARgeted Sequencing](https://github.com/sid-sethi/APTARS): This pipeline processes all the raw sequencing data, does the alignment to a reference genome and collapses reads into unique transcripts. It also annotates transcripts with information such as open reading frame (ORF) predictions and characterises transcript by comparing if they are full splice matches to transcripts in annotation.

[PSQAN (Post Sqanti QC Analysis](https://github.com/sid-sethi/PSQAN): This pipeline applies a set of QC and filtering criteria to remove potential genomic contamination and rare PCR artifacts. Using [`SQANTI3`](https://github.com/ConesaLab/SQANTI3) output of ORF prediction, NMD prediction and structural categorisation based on comparison with the reference annotation, we grouped the identified isoforms into the following categories:

1.  **Non-coding novel** – if predicted to be non-coding and not a full-splice match with the reference
2.  **Non-coding known** – if predicted to be non-coding and a full-splice match with the reference
3.  **NMD novel** – if predicted to be coding & NMD, and not a full-splice match with the reference
4.  **NMD known** – if predicted to be coding & NMD, and a full-splice match with the reference
5.  **Coding novel** – if predicted to be coding & not NMD, and not a full-splice match with the reference
6.  **Coding known (complete match)** – if predicted to be coding & not NMD, and a full-splice & UTR match with the reference
7.  **Coding known (alternate 3’/5’ end)** – if predicted to be coding & not NMD, and a full-splice match with the reference but with an alternate 3’ end, 5’ end or both 3’ and 5’ end.

## Code contents

Within this repository you will find:

<table>
<colgroup>
<col width="11%" />
<col width="88%" />
</colgroup>
<thead>
<tr class="header">
<th>Directory</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><a href="docs" class="uri">docs</a></td>
<td>Contains all <code>.Rmd</code>s and their corresponding <code>.html</code>s describing analyses performed for this project.</td>
</tr>
<tr class="even">
<td><a href="logs" class="uri">logs</a></td>
<td>For any scripts that were run outside of an <code>.Rmd</code> (e.g. scripts from the <a href="scripts" class="uri">scripts</a> directory), a log file was recorded and can be accessed here.</td>
</tr>
<tr class="odd">
<td><a href="raw_data" class="uri">raw_data</a></td>
<td>Data used for the analysis. Most will not be available due to size.</td>
</tr>
<tr class="even">
<td><a href="results" class="uri">results</a></td>
<td>Results from all analyses.</td>
</tr>
<tr class="odd">
<td><a href="scripts" class="uri">scripts</a></td>
<td>Contains analysis scripts. Each script contains a one-line description and is also referenced in its corresponding <code>.Rmd</code>.</td>
</tr>
</tbody>
</table>
