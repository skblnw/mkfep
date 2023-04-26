import os

directories = os.listdir()
both_list = []
free_list = []
complex_list = []

for directory in directories:
    count = 0
    for root, dirs, files in os.walk(directory):
        if "md.xvg" in files:
            count += 1
    if count == 0:
        both_list.append(directory)
    elif count == 24:
        free_count = 0
        complex_count = 0
        for root, dirs, files in os.walk(directory):
            if "md.xvg" in files and "complex" in root:
                complex_count += 1
            elif "md.xvg" in files and "free" in root:
                free_count += 1
        if free_count == 24:
            free_list.append(directory)
        if complex_count == 24:
            complex_list.append(directory)

print("Both:", *both_list)
print("Free:", *free_list)
print("Complex:", *complex_list)
