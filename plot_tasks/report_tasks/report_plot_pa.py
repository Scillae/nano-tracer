from utils.report_plot import pa_3d_report_plot
from utils.tools import chkdir
import os.path

def report_plot_pa():
    conc = 0.5
    temp = 20
    arm_num = 4
    dims_ls = None
    conf_suffix = '-jun_0'

    conc_list = [conc]
    temp_list = [temp]
    arm_num_list = [arm_num]
    data = (conf_suffix, dims_ls, conc_list, temp_list, arm_num_list)
    plt = pa_3d_report_plot(data)
    chkdir(os.path.dirname('report/pa3d.png'))
    plt.savefig('report/pa3d.png',dpi=500)
    plt.close()
    return True