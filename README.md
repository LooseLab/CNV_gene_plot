To initiate a python environment:

```python -m venv cnv_env```

Install package requirements:

```pip install -r requirements.txt```

Activate environment:

```source cnv_venv/bin/activate```


Usage: ```cnv_gene_plot.py [OPTIONS] path/to/BAM```


Options:

  ```-c, --chromosome TEXT```       Chromosome to plot.
  
  ```-g, --gene_panel_path path/to/bed```  Gene bed file to choose.
  
  ```-s, --sub_gene_path path/to/gene_file```    Use file to select genes to display from bed file (.txt format in format of select_genes.txt) 
  
  ```-x, --x_coords min-max```         Hypen-separated base coordinates to plot forthe
                              selected chromosome.
                              
  ```--help```                      Show this message and exit.

Examples:
```
# Generates a CNV plot of chromosome 1 from position 10000000 to 20000000

python3 cnv_gene_plot.py /home/bam_files/test.bam -c chr4 -x 10000000-20000000

# Generates plot of chromosome 1 with all genes from the rCNS panel labelled accordingly

python3 cnv_gene_plot.py /home/thomas/sort_ds1305_Intraop0002_2.htom.bam -c chr4 -g rCNS2_panel_name_uniq.bed 

# Generates a plot of chromosome 1 with selected genes found in select_genes.txt from the rCNS2 panel

python3 cnv_gene_plot.py /home/thomas/sort_ds1305_Intraop0002_2.htom.bam -c chr1 -g rCNS2_panel_name_uniq.bed -s select_genes.txt

```


