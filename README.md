# Sparks Renderer

## Advanced Computer Graphics: Final Project

### Description

The final project is about image creation. You must put
together a coherent scene using geometry, material, and texture that
you made and/or found online. Code your own path tracing algorithm
with acceleration structure to render out a beautiful image.

### TA grading period

Dec 26 — Jan 8

### Grading

Total score: 55 points
- Base: a path tracing algorithm that could handle diffusive material and specular material correctly
with a proper acceleration structure (in total 30 pts)
  - We will provide several standard test cases (~3) to verify the correctness of your algorithm
- Additional features (choose your own but the maximum score won’t exceed 55 pts):
  - Build your own scene: with aesthetics consideration, with geometry you make from scratch and/or
find online with the source highlighted (up to 8 pts)
  - Material: create a (non-trivial) customized material (up to 5 pts)
  - Texture: create your own (non-trivial) texture with the proper texture mapping (up to 5 pts)
  -  Anti-aliasing: implement an anti-aliasing algorithm (up to 2 pts)
  - Importance sampling: using better sampling algorithms for path tracing (up to 5 pts)
  - Simulation based content creation (up to 5 pts)
  - Special material rendering: participating media, hair, translucent material, etc. (up to 5 pts for each)
  - Special visual effect: motion blur, depth of field, etc. (up to 3 pts for each)


## Build The Project

You need to **fork** this repo, and do changes on your forked repo!

You need to **fork** this repo, and do changes on your forked repo!

You need to **fork** this repo, and do changes on your forked repo!

Clone the repo (With submodules recursively):

```
git clone link-to-your-repo.git --recursive
```

or

```
git clone link-to-your-repo.git
cd sparkium
git submodule update --init --recursive
```

### Using Visual Studio 2022

Make sure that you installed Visual Studio with component `Desktop development with C++`.
If not, please open Visual Studio Installer and `Modify` your installation.

The repo contains configuration file `CMakeSettings.json`.
You only need to open the repo folder with Visual Studio on your local machine,
then Visual Studio should start CMake configuration process automatically.

When CMake configuration is successfully finished,
you should be able to select build type and debug target.
You should select `sparks.exe` as the debug target and run the program.

### Using Visual Studio Code

On Windows, you also should install Visual Studio with component `Desktop development with C++`.
We do this for the MSVC compiler, which has the best compatibility on Windows.
MinGW and MSYS2 are not recommended!!!

The repo also contained a vscode profile.
You can just open the repo folder on you local machine and configure CMake through settings on the bottom bar.
Your VSCode should have installed extensions `CMake Tools` and `C/C++` from Microsoft.

To run your program, you should select target `sparks` first. Then, you're recommended to use `Run->Start Debugging` and `Run->Run Without Debugging`, start from the menu bar.

## Project Framework

All the corresponding codes of this project is under `src/sparks`.
You should explore the code and get familiar with the framework.

The following stuff is what you need to know about this framework.

### Scene

The most important part of a renderer is its scene.
A scene contains different kinds of assets.

Codes for describing the scene is under `src/sparks/assets`.

#### Textures

#### Model

##### Mesh

#### Material

#### Entity

### Renderer

#### Path Tracer
