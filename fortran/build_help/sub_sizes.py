import subprocess
import os
from string import Template
import sys


def main():
    print("Running sub_sizes")

    try:
        result = subprocess.check_output(["./sizes"])
    except OSError:
        result = subprocess.check_output(["sizes.exe"])

    split_result = result.split()

    if sys.version_info >= (3, ):
        if split_result[0].strip().decode() != "SIZES":
            raise ValueError("First line is not SIZES")

        subs = {"FOMP_SCHED_KIND": split_result[1].strip().decode(),
                "FOMP_LOCK_KIND": split_result[2].strip().decode(),
                "FOMP_NEST_LOCK_KIND": split_result[3].strip().decode(),
                "FOMP_SCHED_STATIC": split_result[4].strip().decode(),
                "FOMP_SCHED_DYNAMIC": split_result[5].strip().decode(),
                "FOMP_SCHED_GUIDED": split_result[6].strip().decode(),
                "FOMP_SCHED_AUTO": split_result[7].strip().decode()
                }
    else:
        if split_result[0].strip() != "SIZES":
            raise ValueError("First line is not SIZES")

        subs = {"FOMP_SCHED_KIND": split_result[1].strip(),
                "FOMP_LOCK_KIND": split_result[2].strip(),
                "FOMP_NEST_LOCK_KIND": split_result[3].strip(),
                "FOMP_SCHED_STATIC": split_result[4].strip(),
                "FOMP_SCHED_DYNAMIC": split_result[5].strip(),
                "FOMP_SCHED_GUIDED": split_result[6].strip(),
                "FOMP_SCHED_AUTO": split_result[7].strip()
                }

    ompgen_temp_path = os.path.join("..", "ompgen.F90.template")
    ompgen_out_path = os.path.join("..", "ompgen.F90")
    if len(sys.argv) == 3:
        ompgen_temp_path = sys.argv[1]
        ompgen_out_path = sys.argv[2]

    with open(ompgen_temp_path, "r") as ompgen_in:
        ompgen_template = Template(ompgen_in.read())

        ompgen_string = ompgen_template.substitute(subs)

    with open(ompgen_out_path, "w") as ompgen_out:
        ompgen_out.write(ompgen_string)

    print("End sub_sizes")


if __name__ == "__main__":
    main()
