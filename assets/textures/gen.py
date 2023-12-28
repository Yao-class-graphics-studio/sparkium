from PIL import Image
import numpy as np
import math

# 定义图像大小
width, height = 256, 256

# 创建一个空白的RGBA图像数组，初始值为全透明
rgba_array = np.zeros((height, width, 4), dtype = np.uint8)
rgba_array[:, :, 3] = 255

gridR, gridC = 16, 16
gradient = np.zeros((gridR + 1, gridC + 1, 2), dtype = np.float16)

def grad(x: float, y: float):
    vec = [0, 0]
    vec[0] = x * 127.1 + y * 311.7
    vec[1] = x * 269.5 + y * 183.3 

    sin0 = math.sin(vec[0]) * 43758.5453123
    sin1 = math.sin(vec[1]) * 43758.5453123
    
    vec[0] = (sin0 - math.floor(sin0)) * 2.0 - 1.0
    vec[1] = (sin1 - math.floor(sin1)) * 2.0 - 1.0
    length = math.sqrt(vec[0] * vec[0] + vec[1] * vec[1])
    vec[0] /= length
    vec[1] /= length

    return vec

def interpolate(t):
    return 6 * t ** 5 - 15 * t ** 4 + 10 * t ** 3

def perlin(x: float, y: float):
    px, py = math.floor(x), math.floor(y)
    
    g00 = gradient[px, py]
    g01 = gradient[px, py + 1]
    g10 = gradient[px + 1, py]
    g11 = gradient[px + 1, py + 1]
    
    d00 = [x - px, y - py]
    d01 = [x - px, y - py - 1]
    d10 = [x - px - 1, y - py]
    d11 = [x - px - 1, y - py - 1]

    dot00 = d00[0] * g00[0] + d00[1] * g00[1]
    dot01 = d01[0] * g01[0] + d01[1] * g01[1]
    dot10 = d10[0] * g10[0] + d10[1] * g10[1]
    dot11 = d11[0] * g11[0] + d11[1] * g11[1]

    t0 = interpolate(y - py)
    t1 = interpolate(x - px)
    n0 = dot10 * (1 - t0) + dot11 * t0
    n1 = dot00 * (1 - t0) + dot01 * t0
    return n1 * (1 - t1) + n0 * t1

for i in range(gridR + 1):
    for j in range(gridC + 1):
        gradient[i, j] = grad(i, j)

for i in range(width):
    for j in range(height):
        gray = perlin(i / width * gridR, j / height * gridC)
        gray = math.floor((gray + 1) * 255 / 2)
        rgba_array[i, j, 0] = gray

for i in range(width):
    for j in range(height):
        rgba_array[i, j, 1] = rgba_array[height - j - 1, i, 0]
        rgba_array[i, j, 2] = rgba_array[j, height - i - 1, 0]

output_path = 'rgba_image.png'
img = Image.fromarray(rgba_array, 'RGBA')
img.show()
img.save(output_path)