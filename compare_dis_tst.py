import numpy as np
import matplotlib.pyplot as plt

def read_nmnm(fname):
    """Read in the nmnm distribution file, trim the trailing zeros."""
    dt = np.dtype([('nmnm','int'), ('prob','float64')])
    nmnm_file = np.loadtxt(fname, dtype=dt)
    nmnm_file = nmnm_file[ nmnm_file["prob"]>0]
    return nmnm_file


#dt = np.dtype([('nmnm','int'), ('nmnm_prob','float64')])
#nmnm_tst = np.loadtxt('Vsmall.dat', dtype=dt)
#for e in nmnm_tst:
#    print e

def main():
    vsmall = read_nmnm("Vsmall.dat")
    dis_tst = read_nmnm("nmnm_dis_tst.dat")

    fig, ax = plt.subplots()
    vsm = ax.plot(vsmall["nmnm"], vsmall["prob"], 'ro-', lw=4, label='vsmall')
    dis = ax.plot(dis_tst["nmnm"], dis_tst["prob"], 'bs-', lw=2, label='dis_tst')

    plt.show()

if __name__ == "__main__":
    main()
