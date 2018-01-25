import subprocess
import os
from string import Template

def main():
    
    print ("Running sub_sizes")
    try:
        result = subprocess.check_output(["./sizes"])
    except WindowsError:
        result = subprocess.check_output(["sizes.exe"])
    
    split_result = result.split()
    
    if split_result[0].strip() != "SIZES":
        raise ValueError("First line is not SIZES")
    
    subs = {"FOMP_SCHED_KIND" : split_result[1].strip(),
            "FOMP_LOCK_KIND" : split_result[2].strip(),
            "FOMP_NEST_LOCK_KIND" : split_result[3].strip(),
            "FOMP_SCHED_STATIC" : split_result[4].strip(),
            "FOMP_SCHED_DYNAMIC" : split_result[5].strip(),
            "FOMP_SCHED_GUIDED" : split_result[6].strip(),
            "FOMP_SCHED_AUTO" : split_result[7].strip()
            }
    
    
    ompgen_temp_path = os.path.join("..", "ompgen.F90.template")
    ompgen_out_path = os.path.join("..", "ompgen.F90")
    
    with open(ompgen_temp_path, "r") as ompgen_in:
        ompgen_template = Template(ompgen_in.read())
        
        ompgen_string = ompgen_template.substitute(subs)
        
    
    with open(ompgen_out_path, "w") as ompgen_out:
        ompgen_out.write(ompgen_string)
        
    
    print ("End sub_sizes")
if __name__ == "__main__":
    main()