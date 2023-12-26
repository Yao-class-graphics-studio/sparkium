import os
import xpinyin

def convert(str):
    if file.endswith(".mtl"):
        return str
    p = xpinyin.Pinyin()
    return p.get_pinyin(str)

files = os.listdir(".")
print(files)
for file in files:
    with open(file, "r", encoding = "utf-8") as f:
        print(file)
        with open("..\\en-objects\\" + convert(file), "w", encoding = "utf-8") as newf:
            for line in f:
                if line == "":
                    continue
                newf.write(line)