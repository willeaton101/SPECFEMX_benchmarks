

def gen_stn_list(no_stns):
    stn_list = []
    for j in range(no_stns):
        stn_list.append(f"X{str(10 * (j + 1) + 1)}")
    return stn_list