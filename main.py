from plot_tasks.summ_tasks import summ_plot_pa, summ_plot_k2, summ_plot_as, summ_plot_pj, summ_plot_kj, summ_plot_rj, summ_plot_js
from plot_tasks.summ_tasks_juns import summ_plot_pa_jun, summ_plot_k2_jun, summ_plot_as_jun, summ_plot_pj_jun, summ_plot_kj_jun, summ_plot_rj_jun, summ_plot_js_jun


# jobs: arbitrary work batches; should be packed into tasks if reusable.
def summ_plot_main():
    conf_suffix = '' # -jun_10
    dims_ls = [20, 2, 7]
    conc_list = [0.05, 0.1, 0.3, 0.5] # 0.05, 0.1, 0.3, 0.5
    temp_list = [20, 23, 27, 30, 40, 50] # 
    arm_num_list=[3, 4, 5, 6] # 
    task_list = ['Mean', 'ST Dev', 'Skewness'] # m1, std, m
    # summ_list = ['-Patch', '-k2']
    # customization of series ~ conc
    color_list = ['#4994FF','#E55050','#FCC555','#7AA77A'] # np.array((0.1, 0.2, 0.5)).reshape((1,3))  '#4994FF','#E55050','#FFF555','#7AA77A'
    marker_list = ['o','v','^','s'] # 'o','v','^','s'
    # if conf_suffix[:5] == '-jun_':
    #     dims_ls[1] = conf_suffix[-1]

    # summ_plot_pa(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    # summ_plot_k2(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    # summ_plot_as(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)    
    # summ_plot_pj(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    # summ_plot_kj(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    # summ_plot_rj(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    # summ_plot_js(conf_suffix, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    return True

def summ_plot_main_jun():
    conf_suffix = '' # -jun_10
    dims_ls = [20,2,7]
    conc_list = [0.1, 0.5]
    temp_list = [20, 30]
    arm_num_list=[3,4,5,6]
    task_list = ['Mean', 'ST Dev', 'Skewness'] # m1, std, m
    # summ_list = ['-Patch', '-k2']
    # customization of series ~ conc
    color_list = ['#4994FF','#E55050'] # np.array((0.1, 0.2, 0.5)).reshape((1,3))
    marker_list = ['o','v']

    jun_list = [0,1,2,5,10]

    # summ_plot_pa_jun(jun_list, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)   
    # summ_plot_k2_jun(jun_list, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    # summ_plot_as_jun(jun_list, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    # summ_plot_pj_jun(jun_list, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    # summ_plot_kj_jun(jun_list, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    # summ_plot_rj_jun(jun_list, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    summ_plot_js_jun(jun_list, dims_ls, conc_list, temp_list, arm_num_list, task_list, color_list, marker_list)
    # TODO: in the nanostar generation, the center is somehow not correctly stored: always 2 bases.
    return True

def misc():
    # ns_pa_plot(single=True,arms=5,temp=23,conc=0.3)
    # for conc in conc_list:
    #     ns_k2_plot(single=True,arms=4,temp=30,conc=conc)
    dims_ls = [20,2,7]
    arms = 3
    temp = 30
    conc = 0.3
    conf_suffix = ''
    sp_suffix = ''
    label = f'{arms}arms@({temp}C,{conc}M){conf_suffix}{sp_suffix}'
    loose_lbl = f'{temp}C-{conc}M-GPU{sp_suffix}'
    p = f'data/composed_traj/{arms}arms{conf_suffix}/{loose_lbl}/{label}.k2tp'
    # finding_extreme_k2_vals(p) # deprecated. saved in DNA3
    top_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}/{arms}arm-rods-clustered{conf_suffix}.top'
    traj_path = f'../../ox-sync/simul-inputs-{arms}arms{conf_suffix}/{loose_lbl}/trajectory.dat'
    from calc_tasks.patch_angle_calc import patch_angle_calc as pa
    pa(top_path,traj_path,arms,dims_ls)
    return True
# jobs end

if __name__ == '__main__':
    summ_plot_main()
    summ_plot_main_jun()
    # misc()
    print('DONE')