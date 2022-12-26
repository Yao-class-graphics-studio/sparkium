# ���������ļ�˵��

## ���

���������ļ���ʾ����`base.xml`����һ��xml�ļ����ļ���ʽ������Ⱦ��[mitsuba](https://mitsuba2.readthedocs.io/en/latest/src/getting_started/file_format.html)��ʹ�õ�����xml�ļ���

## �ļ���ʽ

��������`scene`�ڻᶨ��������Ԫ�أ����`camera`���ƹ�`envmap`������`mesh`������͵ƹ���һ�������б���ģ����Դ�������Ĭ�ϲ�������ʹ�����ļ���û�ж���Ҳ����Ĭ��ֵ��

�����ļ��Ļ�����ʽΪ��
```
<type name="xxx" value="xxx"/>
```
����typeΪ�������ͣ�nameΪ��������valueΪֵ

���磺

```
<float name="radius" value="0.5"/>
```

��ʾfloat����radius��ֵΪ0.5��

### ���`camera`

| ������       | ����      | ����                      |
| ------------ | --------- | ------------------------- |
| -            | transform | �����λ�ˣ����transform |
| fov          | float     | δʵ��                    |
| focal_length | float     | ���࣬δʵ��              |

### �ƹ�`envmap`

Ԥ���廷��������Ϊ`type=hdr`����Ҫ����ʵ�ֱ��

| ������ | ����    | ����                  |
| ------ | ------- | --------------------- |
| -      | texture | �����⣬ͨ����hdrͼƬ |

###����`mesh`

mesh���Զ���`.obj`��ʽ�ļ���Ҳ����ʹ�ò�������Ԥ��������`Sphere`����`Cube`���ù���ͨ���༭mesh����������`type`ȷ�������ṩ��ͬ�Ĳ��������磺

```
<mesh type="obj">
	<string name="filename" value="../../meshes/plane.obj"/>
</mesh>

<mesh type="Sphere">
	<vec3 name="center" value="0.0 0.0 0.0"/>
	<float name="radius" value="0.5"/>
</mesh>
```

Ϊ���ֲ�ͬ���͵�������롣���±�ʾ��ͬtype��mesh����Ĳ�ͬ����������

| type   | ������   | ����   | ����                     |
| ------ | -------- | ------ | ------------------------ |
| obj    | filename | string | .obj�ļ���·��           |
| Sphere | center   | vec3   | �������λ��             |
|        | radius   | float  | ��İ뾶                 |
| Cube   | center   | vec3   | �����������λ�ã�δʵ�� |
|        | size     | float  | ������Ĵ�С��δʵ��     |

�����Ƕ����������������ͳһ�ı�����

| ������ | ����      | ����                     |
| ------ | --------- | ------------------------ |
| -      | transform | ���嵼�����Ҫ���еı任 |
| -      | material  | ����Ĳ��ʣ����material |

### �任`transform`

�任����ͨ����һ����������������4x4����ͨ��������ת��ƽ�������֡������ǲ�ͬtype��transform���費ͬ�Ĳ���������

| type       | ������      | ����  | ����                   |
| ---------- | ----------- | ----- | ---------------------- |
| lookat     | eye         | vec3  |                        |
|            | center      | vec3  |                        |
|            | up          | vec3  |                        |
| translate  | translation | vec3  | ƽ������               |
| axis-angle | axis        | vec3  | ������ת���ᣬδʵ��   |
|            | angle       | float | ������ת�ĽǶȣ�δʵ�� |
| mat4       | matrix      | mat4  | 4x4����������δʵ��  |

### ����`material`

���ʵ����Կ���ͨ������material��`type`���ԣ���ѡ����`lambertian`��`specular`��`transmissive`��`principled`���ֱ��Ӧ������������������ԡ�

| ������       | ����               | ����                                                   |
| ------------ | ------------------ | ------------------------------------------------------ |
| albedo       | texture��rgb��rgba | ���������albedo�������ǵ�ɫrgb/rgba��Ҳ������������ͼ |
| albedo_color | vec3               | ������ɫ��δʵ��                                       |

###����`texture`

Ŀǰֻʵ���˵���ͼƬ��checkerֻ��Ϊ�˾�����

| type    | ������   | ����   | ����                       |
| ------- | -------- | ------ | -------------------------- |
| bitmap  | filename | string | ͼƬ·��                   |
| checker | color1   | vec3   | ���̸��һ����ɫ��δʵ��   |
|         | color2   | vec3   | ���̸����һ����ɫ��δʵ�� |



`base.xml`���Ը������������ļ���

```
<?xml version="1.0" encoding="utf-8"?>
<scene>
	
	<camera>
		<transform type="lookat">
			<vec3 name="eye" value="2.0 1.0 3.0"/>
			<vec3 name="center" value="0.0 0.0 0.0"/>
			<vec3 name="up" value="0.0 1.0 0.0"/>
			<!--<lookat eye="2.0 1.0 3.0" center="0.0 0.0 0.0" up="0.0 1.0 0.0"/>-->
		</transform>
	</camera>
	
	<envmap type="hdr">
		<string name="filename" value="../../textures/envmap_clouds_4k.hdr"/>
	</envmap>
	
	<mesh type="obj">
		<string name="filename" value="../../meshes/plane.obj"/>
		
		<material type="lambertian">
			<rgb name="albedo" value="0.8 0.8 0.8"/>
		</material>
	</mesh>

	<mesh type="Sphere">
		<vec3 name="center" value="0.0 0.0 0.0"/>
		<float name="radius" value="0.5"/>
		
		<transform type="translate">
			<vec3 name="translation" value="0.0 0.5 0.0"/>
		</transform>
		
		<material type="lambertian">
			<texture name="albedo">
				<string name="filename" value="../../textures/earth.jpg"/>
			</texture>
		</material>
	</mesh>
	
</scene>
```

