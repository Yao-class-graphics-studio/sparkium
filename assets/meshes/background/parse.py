import os

output = []
path_pre = "../../" 

def add_tab(arr): #缩进
    ans = arr
    for i,_ in enumerate(ans):
        ans[i] = "\t"+ans[i]
    return ans

def parse_mtl(filename): #texture 存放在 textures/background/ 子目录下
    ans=[]
    tmp=[]
    have_diffuse = True
    flag = False
    with open(filename, "r", encoding="utf-8") as f:
        for st in f:
            if st == "":
                continue
            content = st.split()
            if(len(content)==0):
                continue
            if content[0] == "newmtl" and content[1] == "金边":
                flag |= True
            if(content[0]=="map_Kd"):
                texture_filename = content[1].split('\\')[-1]
                tmp.append('<albedo_texture value="'+path_pre+'textures/background/'+texture_filename+'"/>')
            elif(content[0]=="Kd"):
                color = ""
                for i in range(3):
                    color += content[i+1]
                    if(i!=2):
                        color += ' '
                tmp.append('<albedo value="'+color+'"/>')
    if flag == True:
        tmp.append('<metallic="1"/>')
        tmp.append('<roughness="1"/>')
        ans.append('\t<material type="principled">')
        ans += add_tab(add_tab(tmp))
    else:
        if not have_diffuse:
            tmp.append('<albedo value="1 1 1"/>')
        ans.append('\t<material type="lambertian">')
        ans += add_tab(add_tab(tmp))
    ans.append('\t</material>')
    return ans

def parse_obj(filename): # obj 也是存放在 meshes/background/ 子目录下
    ans = []
    ans.append(r"""<model type="obj">""")
    ans.append('\t<filename value="'+path_pre+"meshes/background/"+filename+'"/>')

    # ans.append('\t<material type="lambertian">')
    with open(filename,"r", encoding = "utf-8") as f:
        while True:
            s = f.readline()
            if(s==""):
                break
            content = s.split()
            if(len(content)==0):
                continue
            if content[0]=="mtllib":
                ans += parse_mtl(content[1])
            else:
                continue
    # ans.append('\t</material>')
    ans.append('</model>')
    return ans

files = os.listdir(".")
for file in files:
    if (not file.endswith(".obj")):
        continue
    print(file)
    output += add_tab(parse_obj(file))
    output.append('\n')

with open("tmp.xml","w", encoding = "utf-8") as f:
    f.write('<?xml version="1.0" encoding="utf-8"?>'+'\n')
    f.write('<scene>'+'\n')
    for s in output:
        f.write(s+'\n')
    f.write('</scene>'+'\n')