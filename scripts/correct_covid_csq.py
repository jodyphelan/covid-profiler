import sys
import re

for l in sys.stdin:
    row = l.rstrip().split("\t")
    if "part2" in l:
        info = row[7].split(";")
        for i in range(len(info)):
            if info[i][:4]=="BCSQ":
                new_bcsq = []
                for bcsq in info[i].split(","):
                    if bcsq[5]=="@":
                        new_bcsq.append(bcsq)
                    else:
                        tmp = bcsq.split("|")

                        re_obj1 = re.search("(\d+)([A-Z]+)>(\d+)([A-Z]+)",tmp[5])
                        re_obj2 = re.search("(\d+)([A-Z]+)",tmp[5])
                        if re_obj1:
                            tmp[5] = "%s%s>%s%s" % ((int(re_obj1.group(1)) + 4401), re_obj1.group(2), (int(re_obj1.group(3)) + 4401), re_obj1.group(4))
                        elif re_obj2:
                            tmp[5] = "%s%s" % ((int(re_obj2.group(1)) + 4401), re_obj2.group(2))
                        else:
                            quit(l)
                        new_bcsq.append("|".join(tmp))
                info[i] = ",".join(new_bcsq)
        row[7] = ";".join(info)
    sys.stdout.write("\t".join(row)+"\n")
