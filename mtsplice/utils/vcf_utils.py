HEADER = """##fileformat=VCFv4.0
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##contig=<ID=MT,length=16569>
"""

COLS = '''#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO'''

# name convertor


def variant_name_convert(var):
    '''convert variant name
    e.g. 10:100007693:T:['C'] --> 10_100007693_T_C or the other way around
    1:897349:G>A --> 1_897349_G_A
    '''
    if "['" in var:
        var_split = var.split(':')
        var_split[-1] = var_split[-1].strip("['").strip("']")
        return '_'.join(var_split)
    elif '_' in var:
        var_split = var.split("_")
        return var_split[0] + ":" + var_split[1] + ":" + var_split[2] + ":" + "['" + var_split[3] + "']"
    elif ">" in var:
        var_split = var.split(":")
        ref_var = var_split[-1].split(">")
        return '_'.join(var_split[:-1] + ref_var)
