file = open("domains_homepage")

for line in file:
    strline = str(line).strip()
    if strline.split()[8] == "FM":
        print strline.split()[5][:-1], strline.split()[7], strline.split()[6]
file.close()
