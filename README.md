This tool now includes ruptures changepoint detection using the ```BinSeg``` algorythm, for estimation of CNV breakpoints. Happy plotting.

The centromeric/acroscentric regions are identified using the t2t-chm13-v1.0 database for GRCh38 found [here.](https://genome.ucsc.edu/cgi-bin/hgTables?db=hub_2395475_t2t-chm13-v1.0&hgta_group=cenTelo&hgta_track=hub_2395475_cenRegions&hgta_table=hub_2395475_cenRegions&hgta_doSchema=describe+table+schema)


## Installation
To initiate a python environment:

```python -m venv cnv_env```

Install package requirements:

```pip install -r requirements.txt```

Activate environment:

```source cnv_venv/bin/activate```

## Usage
```
Usage: cnv_gene_plot_ruptures.py [OPTIONS] path/to/BAM


Options:

 -c, --chromosome TEXT     Chromosome to plot.
 
 -g, --gene_panel_path path/to/bed.bed  Gene bed file to choose
 
 -s, --sub_gene_path path/to/genes.txt   Use file to select genes to display
 
 -x, --x_coords INTEGER-INTEGER     Hyphen-separated base coordinates to plot for
                              the selected chromosome
                              
 -o, --output_file TEXT.pdf   Select output PDF file
 
 -r, --ruptures             Include ruptures changepoint detection
 
  --help                    Show this message and exit.
```

## Examples
```sh
# Generates a CNV plot of chromosome 1 from position 10000000 to 20000000

python3 cnv_gene_plot.py /home/bam_files/test.bam -c chr4 -x 10000000-20000000

# Generates plot of chromosome 1 with all genes from the rCNS panel labelled accordingly

python3 cnv_gene_plot.py /home/thomas/sort_ds1305_Intraop0002_2.htom.bam -c chr4 -g rCNS2_panel_name_uniq.bed 

# Generates a plot of chromosome 1 with selected genes found in select_genes.txt from the rCNS2 panel with ruptures changepoint detection

python3 cnv_gene_plot.py /home/thomas/sort_ds1305_Intraop0002_2.htom.bam -c chr1 -g rCNS2_panel_name_uniq.bed -s select_genes.txt -r 

```


