import pandas
import numpy
import matplotlib.pyplot as plt

def load_expr(array_exp_file, num_sel, delimiter='\t'):
    """Load the dataset in a pandas dataframe"""

    if os.path.isfile(array_exp_file):
        expr = pandas.read_csv(array_exp_file, sep=delimiter)
    (rows, cols) = expr.shape
    if num_sel > cols:
        num_sel = cols
    expr_sel = expr.ix[:,0:num_sel]
    expr_test = expr.ix[:,num_sel::]
    return expr_sel

def get_coef_var(expr, n):
    """
    Determine the coefficient of variation for each gene
    and return the top N genes
    """

    cv = expr.mean(axis=1)/expr.std(axis=1)
    cv.sort()
    cv_sel = cv[:n]
    return cv_sel


def get_updated_expr(expr, cv_sel):
    """ build a new matrix with selected genes"""

    rows, cols = expr.shape
    for i in xrange(rows):
        if i not in cv_sel.keys():
            expr = expr.drop(i, axis=0)
    return expr

#create a heatmap of that matrix

def create_heatmap(expr):
    plt.pcolor(expr)
    plt.savefig("heatmap.png")
    plt.clf()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ccle.py', description='Generate heatmaps')
    parser.add_argument("--expr", type=str, help='Path to expression dataset', required=True)
    parser.add_argument("--n_samples", type=str, help='No of samples selected', required=True)
    parser.add_argument("--n_features", type=str, help='No of features selected', required=True)

    args = parser.parse_args()

    expr = load_expr(args.expr, args.n_samples)
    cv_sel = get_coef_var(expr, args.n_features)
    expr = get_updated_expr(expr, cv_sel)
    create_heatmap(expr)


