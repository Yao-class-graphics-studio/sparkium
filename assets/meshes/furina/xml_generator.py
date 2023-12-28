import os

output = []
path_pre = "../../"

def add_tab(arr):
    ans = arr
    for i,_ in enumerate(ans):
        ans[i] = "\t"+ans[i]
    return ans

def parse_mtl(filename):
    ans=[]
    have_diffuse = True
    with open(filename, "r", encoding="utf-8") as f:
        while True:
            st = f.readline()
            if(st==""):
                break
            content = st.split()
            if(len(content)==0):
                continue
            if(content[0]=="map_Kd"):
                texture_filename = content[1].split('\\')[-1]
                ans.append('<albedo_texture value="'+path_pre+'textures/furina/'+texture_filename+'"/>')
            elif(content[0]=="Kd"):
                color = ""
                for i in range(3):
                    color += content[i+1]
                    if(i!=2):
                        color += ' '
                ans.append('<albedo value="'+color+'"/>')
    if not have_diffuse:
        ans.append('<albedo value="1 1 1"/>')
    return ans

def parse_obj(filename):
    ans = []
    ans.append(r"""<model type="obj">""")
    ans.append('\t<filename value="'+path_pre+"meshes/furina/"+filename+'"/>')

    ans.append('\t<material type="lambertian">')
    with open(filename,"r", encoding = "utf-8") as f:
        while True:
            s = f.readline()
            if(s==""):
                break
            content = s.split()
            if(len(content)==0):
                continue
            if content[0]=="mtllib":
                ans += add_tab(add_tab(parse_mtl(content[1])))
            else:
                continue
    ans.append('\t</material>')
    ans.append('</model>')
    return ans

files = os.listdir(".")
for file in files:
    if (not file.endswith(".obj")):
        continue
    output += add_tab(parse_obj(file))

with open("tmp.xml","w", encoding = "utf-8") as f:
    for s in output:
        f.write(s+'\n')
