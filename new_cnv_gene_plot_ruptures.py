import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import click
from cnv_from_bam import iterate_bam_file
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
import ruptures as rpt

@click.command()
@click.argument('input_dir', type=click.Path(exists=True))
@click.option('-c', '--chromosome', default='chr1', help='Chromosome to plot.')
@click.option('-g', '--gene_panel_path', default='', help='Gene bed file to choose')
@click.option('-s', '--sub_gene_path', help='Use file to select genes to display')
@click.option('-x', '--x_coords', default='', help='Hyphen-separated base coordinates to plot for the selected chromosome')
@click.option('-o', '--output_file', default='output.pdf', help='Select output PDF file')
@click.option('-r', '--ruptures', is_flag=True, show_default=True, default=False, help='Include ruptures changepoint')
def read_bam(input_dir, chromosome, gene_panel_path, sub_gene_path, x_coords, output_file, ruptures):
    bam_files = [f for f in os.listdir(input_dir) if f.endswith('.bam')]
    bam_files = [os.path.join(input_dir, f) for f in bam_files]

    if not bam_files:
        print(f"No BAM files found in directory: {input_dir}")
        return

    with PdfPages(output_file) as pdf:
        for file_path in bam_files:
            bam_path = Path(file_path)
            
            bam_name_suff = os.path.basename(file_path)
            bam_name = bam_name_suff.split('.')[0]

            if chromosome == 'chrx':
                chromosome = 'chrX'
            if chromosome == 'chry':
                chromosome = 'chrY'

            # Parse the x_coords option
            x_start, x_end = None, None
            if x_coords:
                try:
                    x_start, x_end = map(int, x_coords.split('-'))
                except ValueError:
                    raise ValueError("Invalid format for x_coords. Use hyphen to separate start and end coordinates, e.g., '1000000-2000000'.")

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 5), dpi=800, constrained_layout=True)
            ax.set_ylim((0, 5))
            ax.xaxis.set_major_locator(plt.MaxNLocator(10))

            sub_gene_list = []
            bed_array = None

            if gene_panel_path:
                try:
                    gene_bed_df = pd.read_csv(gene_panel_path, sep='\t', header=None)
                    bed_array = gene_bed_df.to_numpy()
                except FileNotFoundError:
                    print(f"File not found: {gene_panel_path}. Continuing without gene panel plot.")

            if sub_gene_path:
                with open(sub_gene_path, 'r') as sub_gene_file:
                    sub_gene_list = [line.strip() for line in sub_gene_file]

            result = iterate_bam_file(bam_path, _threads=4, mapq_filter=60)
            bin_width = result.bin_width

            for contig, cnv in result.cnv.items():
                if contig == chromosome:
                    x_values = np.arange(len(cnv)) * bin_width
                    cnv = np.array(cnv)

                    if x_start is not None and x_end is not None:
                        mask = (x_values >= x_start) & (x_values <= x_end)
                        x_values = x_values[mask]
                        cnv = cnv[mask]

                    ax.scatter(x=x_values, y=cnv, s=4)

                    if ruptures:
                        # Perform changepoint detection
                        algo = rpt.Binseg(model="l2").fit(cnv)
                        result = algo.predict(n_bkps=5)

                        # Display changepoints
                        for bkpt in result:
                            if bkpt < len(x_values):
                                ax.axvline(x=x_values[bkpt], color='green', linestyle='--')

            if bed_array is not None:
                for array in bed_array:
                    if sub_gene_path:
                        if array[0] == chromosome:
                            if array[3] in sub_gene_list:
                                start = array[1]
                                end = array[2]
                                label = str(array[3])
                                if (x_start is None and x_end is None) or (start >= x_start and end <= x_end):
                                    ax.axvspan(start, end, alpha=0.5, color='red')
                                    ax.text((start + end) / 2, 5.1, label, rotation=45, verticalalignment='bottom', fontsize=8, color='red')

                    else:
                        if array[0] == chromosome:
                            start = array[1]
                            end = array[2]
                            label = str(array[3])
                            if (x_start is None and x_end is None) or (start >= x_start and end <= x_end):
                                ax.axvspan(start, end, alpha=0.5, color='red')
                                ax.text((start + end) / 2, 5, label, rotation=45, verticalalignment='bottom', fontsize=8, color='red')

            title = f"{bam_name}: {chromosome}"
            if x_coords:
                title += f" {x_coords}"
            plt.title(title, pad=80)

            pdf.savefig(fig)
            plt.close(fig)

if __name__ == '__main__':
    read_bam()

