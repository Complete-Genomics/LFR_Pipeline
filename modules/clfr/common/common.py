"""
commonly used functions TODO: combine to utility.py
"""
import sys
from PIL import Image
import math

def combine_plots(frag_len_min, frag_len_max, bin_frag_len, dirname, type, library_type):
    """
    combine mean_cov_bin_boxplot into 1 plot, input results .png from frag_mean_cov_boxplot.v2.py
    """
    images = []
    try:
        for x in range(frag_len_min, frag_len_max+bin_frag_len, bin_frag_len):
            images.append(Image.open(dirname+'/{}/mean_cov_bin_boxplot_'.format(type)+str(x)+'.png'))
    except Exception as e: print(e) 
    widths, heights = zip(*(i.size for i in images))
    ## 0.3-3k(28) n_row, n_col=4*7, 1-5.1k(42) =6*7
    num_sample = len(images)
    n_row = int(num_sample//math.sqrt(num_sample))
    if num_sample%n_row > 0:
        n_col = num_sample//n_row +1
    else:
        n_col = num_sample//n_row
    # if library_type =='gdna':
    #     n_row, n_col = 6,7 
    # elif library_type =='mrna':
    #     n_row, n_col = 4,7

    total_width = max(widths)*n_row
    max_height = max(heights)*n_row

    new_im = Image.new('RGB', (total_width, max_height), color="white")
    y_offset = images[0].size[1]
    x_offset = images[0].size[0]
    n = int((1+(frag_len_max-frag_len_min)/bin_frag_len)/n_col)+1
    j=0   
    while j<n_row:
        for i in range(j*n_col,(j+1)*n_col):
            if i<len(images):
                new_im.paste(images[i], (x_offset*(i-j*n_col),y_offset*j))
        j+=1    

    new_im.save(dirname+'/{}/mean_cov_bin_boxplot_combined.jpg'.format(type))

def bed_sum(bed_file_merged):
    total_len_merged = 0
    with open(bed_file_merged, 'r') as f:
        for line in f:
            info = line.strip().split('\t')
            chrom = info[0]
            start = int(info[1])
            end = int(info[2])
            total_len_merged +=end-start+1
    print(total_len_merged)
    return total_len_merged

if __name__ == "__main__":
    frag_len_min = 300
    frag_len_max = 3000
    bin_frag_len = 100
    library_type = 'gdna'
    dirname= 'Calc_Frag_Length_50000/'
    type='bin_num10'
    combine_plots(frag_len_min, frag_len_max, bin_frag_len, dirname, type, library_type)
    type='bin_size100'
    combine_plots(frag_len_min, frag_len_max, bin_frag_len, dirname, type, library_type)