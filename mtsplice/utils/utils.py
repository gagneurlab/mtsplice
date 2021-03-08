import os, sys
import matplotlib.pyplot as plt
import numpy as np

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
        
def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')



def row_mean_impute(X, axis=-1):
    return np.where(np.isnan(X), np.nanmean(X, axis=axis, keepdims=True), X)
    

def get_var_side(var):
    ''' Get exon variant side '''
    varstart, ref, alt, start, end, strand = var

    # for long insertion or deletion which can go out side of exon
    # CA->CAGG
    varend = varstart + max(len(ref), len(alt)) - 1

    # left normalization
    for i in range(min(len(ref), len(alt))):
        if ref[i] == alt[i]:
            varstart += 1
        else:
            break

    if strand == "+":
        if varstart < start-50:
            return "left intron"
        elif varend <= start+3 and varstart >= start-50:
            return "acceptor"
        elif varend > end+13:
            return "right intron"
        elif varend <= end+13 and varstart >= end-5:
            return "donor"
        else:
            return "exon"
    else:
        if varstart < start-13:
            return "right intron"
        elif varstart >= start-13 and varend <= start+5:
            return "donor"
        elif varend > end+50:
            return "left intron"
        elif varstart >= end-3 and varend <= end+50:
            return "acceptor"
        else:
            return "exon"    


def get_cor_region(var):
    '''Get variants in the core region
    exon + donor and acceptor covered region
    '''
    varstart, ref, alt, start, end, strand = var
    
