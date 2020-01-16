The formats and utilities documented here were developed at the [Carnegie Mellon Textiles Lab](https://textiles-lab.github.io/) as part of the [Visual Knit Programming](http://visual.knit.zone) project and further enhanced for use in [KnitDB](http://db.knit.zone) and ongoing research.
For bugs, file a github issue.
For other inquires, please contact [Jim McCann](http://www.cs.cmu.edu/~jmccann/).

# `.smobj`: a format for augmented-stitch-mesh-like things

Augmented Stitch Meshes represent knit structures and their dependencies by embedding yarns in the faces of a (tri/quad/pentagon/etc) mesh, and associating types and directions with each edge.

The `.smobj` format stores an augmented stitch mesh as follows:

```
#text format, like .obj with...
#library of face types as name followed by edge types (+/- indicate [optional] edge direction; edge labels are arbitrary strings)
L knit-to-right -l +y +l -y
L knit-to-right -l2 +y +l2 -y #face names include the edge types, so this is different than the previous face
#vertices as X Y Z:
v 1.0 2.2 1.0
v 1.0 1.0 1.0
v 3.0 2.2 1.0
v 3.0 1.0 1.0
#faces as 1-based vertex index lists (and, possibly, texture coordinates):
f 2 4 3 1
#and, for each face, a type from the library (1-based index into library):
T 1
#(optional) 1-based source line number of knitout file that made this face (use 0 if line not known):
N 15
#list of connections between faces:
# face#/edge#, both one-based; negative edges imply *reversing* the edge
# this allows non-orientable connections, needed for knitout
# (NOTE: on a nicely oriented mesh, all 'e' commands will have one negative edge)
e 1/1 2/-4
#(optional) list of checkpoints on yarns and desired length between them:
# checkpoint list starts with unit library:
U 1 1.0 #unit definition -- there is a length unit called '1' with [default] length '1.0'
#By convention, the first unit will always be called '1' and will always have length '1.0'
U n 0.4 #unit definition -- there is a length unit called 'n' with [default] length '0.4'
U s65 0.2 #unit definition -- there is a length unit called 's65' with [default] length '0.2'
#checkpoints are listed *in order* along yarns
#  any yarn with checkpoints *MUST* be covered from start-to-end
c 1/1/1 0.1 1 1.0 2 2.0 3 #checkpoint at face/edge/yarn crossing with 0.1 + 1.0*n + 2.0*s65 length following
c 1/3/1 #last checkpoint on a yarn will always have zero following length
#NOTE: checkpoints define yarn orientation; yarns without checkpoints will be arbitrarially oriented by smobj-to-yarns
```

# `.yarns`: a format for yarn curves

The yarn format stores one or more yarns as polylines.
Typically, these polylines are interpreted for viewing by connecting subsequent midpoints with quadratic bezier curves in order to create a smooth path.

A `.yarns` file consists of several chunks, where each chunk contains a 4-byte header, a 4-byte size (length of chunk data in bytes), and a flat array of uniform data elements.

In other words, each chunk corresponds to a `vector< T >` for some `T`.

The specific chunks and their order:
```
(1) Points
  Header: 'f3..'
  Size: 12*
  Contents: N point locations stored as floating point numbers:
    struct {
	  float x;
	  float y;
	  float z;
	}

(2) Source Line Numbers:
  Header: 'src.'
  Size: 4*N
  Contents: N source line numbers (same length as points) stored as 32-bit unsigned integers
  NOTE: line numbers correspond to the line number for the smobj face that created the next segment. Line numbers for the last point in a yarn are set to the same value as the second-to-last by convention.

(3) Yarn Information
  Header: 'yarn'
  Size: 16*N
  Contents: N yarn descriptors stored as:
    struct {
	  uint32_t point_begin; //index of first point in yarn
	  uint32_t point_end; //index of last+1 point in yarn
	  //Note that yarns will *not* overlap in the points list
	  float radius; //radius of yarn
	  uint8_t r,g,b,a; //color for yarn
	}
[remaining sections are work-in-progress]
(4) Strings
  Header: 'strs'
  Size: N
  Contents: N bytes of string data
(5) Units
  Header: 'unit'
  Size: 12*N
  Contents: N unit descriptors stored as:
    struct {
      uint32_t name_begin; //index of first character of unit name
      uint32_t name_end; //index of last character of unit name
      float length; //default/reccomended length of given unit <-- might remove this
    }
  NOTE: the first unit will always be named '1' and will always represent absolute lengths.
(6) Checkpoints:
  Header: 'chk.'
  Size: 12*N
  Contents: N checkpoints stored as:
    struct {
      uint32_t point; //point to add slack after
      float length; //amount of length to add
      uint32_t unit; //unit of length
    }
  NOTE: checkpoints must be stored sorted by (point, unit)
```

Note that each "Checkpoint" structure adds yarn length between that checkpoint and the next *checkpoint* (or the end of the yarn) -- this is different than adding length between checkpoints and their next yarn point.
This ambiguity is required because the allocation of length to yarn segments is generally not known.


# Utilities

This repository contains a few utilities that make it easy to work with .smobj files.

Build on OSX or Linux using `make` in the `utilities` directory.
You will need glm installed in a system-wide path (or to modify the Makefile as is proper).

Build on Windows using `nmake -f Makefile.windows` from the Visual Studio Command Prompt in the `utilities` directory.
You will need a checkout of glm in the `utilities\glm` directory.
Windows builds are tested less frequently than OSX or Linux builds so you may need to fix a few warnings.

## .knitout to .smobj (Work in Progress)

The `knitout-to-smobj` utility converts knitout instructions to an smobj description.
```
./knitout-to-smobj <input.knitout> <output.smobj>
```

Given the nature of the conversion, the output file often contains elongated yarns, though the yarn length checkpoints it produces should provide a more reasonable notion of how much they need to be shrunk.

### TODO
 - Testing required, especially of distance checkpoints during plating.
 - Yarn bring-in faces should probably be hoisted.
 - Handling of distance checkpoints around splits is not ideal.

## .smobj to .yarns (Work in Progress)

The `smobj-to-yarns` utility loads a face library and smobj file and exports a `.yarns` file which describes the yarn paths described by the file.
```
./smobj-to-yarns <input.smobj> <input-library.sf> <output.yarns>
```

It uses the notion of generalized barycentric coordinates to warp the face from the library.

### TODO
 - Testing required, especially of new features.

# Viewers

## web/view-yarns.html (Work in Progress)

Viewer for yarns files built with javascript and WebGL.
To use, open the file in a browser and drag in (or click to open) a ```.yarns``` file.

### TODO
 - Units panel should update unit lengths.
 - Toggle viewing mode between yarn color and stretch color.
 - Some way to read back line numbers (mouse-over?)


# Standard Face and Edge Types (Knitout)

This section describes the face library our knitout-to-smobj code uses to represent the result of machine knitting.
This library is stored in the `faces/knitout.sf` file.

The library includes `yN` (N yarn), `lN` (N loop), and 'x' (empty) edges.

Here are the face types; edges are named in CCW order from the bottom left:
```
#Basic machine operations:
#front/back knit:
L knit-to-left -l1 -y1 +l1 +y1
L knit-to-right -l1 +y1 +l1 -y1
L purl-to-left -l1 -y1 +l1 +y1
L purl-to-right -l1 +y1 +l1 -y1
#(and versions that take lN -> lM via yM)

#Tuck
L tuck-behind-to-left -l0 -y1 +l1 +y1
L tuck-behind-to-right -l0 +y1 +l1 -y1
L tuck-infront-to-left -l0 -y1 +l1 +y1
L tuck-infront-to-right -l0 +y1 +l1 -y1
#(and versions that take lN -> l(N+M) via yM)

#Split
# top and bottom edge have two loop segments, connected to needle pair being split between
# by convention, the left segment is the front loop and the right is the back loop
L split-front-to-right -l1 -l1 +y1 +l1 +l2 -y1
L split-front-to-left -l1 -l1 -y1 +l1 +l2 +y1
#(and versions that take -lX -lY to +lZ +l(X+Y) via yZ)
L split-back-to-right -l1 -l1 +y1 +l2 +l1 -y1
L split-back-to-left -l1 -l1 -y1 +l2 +l1 +y1
#(and versions that take -lX -lY to +l(X+Y) +lZ via yZ)

#loop / yarn routing:
L loop -l1 x +l1 x
L yarn-to-right x +y1 x -y1
L yarn-to-left x -y1 x +y1
#yarn enters going left/right into diagonal bottom face, leaves going left/right from diagonal top face:
L yarn-left-up-right -y1 x +y1 x
L yarn-left-up-left -y1 x +y1 x
L yarn-right-up-right -y1 x +y1 x
L yarn-right-up-left -y1 x +y1 x
#yarn plating (lower edge connects to bottom of stack so stack [1 2] splits to lower edge [1] and upper edge [2]):
L yarn-plate-to-right x +y2 x -y1 -y1
L yarn-plate-to-left x -y1 -y1 x +y2
#also +y(N+M) from -yN, -yM
L yarn-unplate-to-right x +y1 +y1 x -y2
L yarn-unplate-to-left x -y2 x +y1 +y1
#also +yN +yM from -y(N+M)

```




# Notes and Work-In-Progress Stuff

## Deprecated Things

Our format used to also contain the following, but they have been replaced:

```
# ---- old data ----
#or (2) one of 'tk', 'te', 'ts' to label edges where 0 = loopwise and 1= yarnwise
tk 1 0 1 0 
# and st to mark short-row faces, unused mostly
st 14
# 's' specifies starting face(s)
s f#
```

## Extension: Texture Coordinates

One may add texture coordinates to smobj files as follows:
```
#(optional) texture coorinates as vt commands:
vt 0.0 0.0
vt 1.0 0.0
vt 1.0 1.0
vt 0.0 1.0
#(optional) extra layers of texture coordinates (same length as first list):
vt2 0.5 0.5
vt2 1.0 0.5
vt2 1.0 1.0
vt2 0.5 1.0
#faces include texture coordinates [as per .obj]:
f 2/1 4/2 3/3 1/4
#(optional) indicate which image corresponds to which texture coordinate:
tex foo.jpg
tex2 bar.jpg
```



## What /Should/ Be In The File

 - A collection of faces.
 - For each face, a type (which also implies a type for all edges)
 - Edge-to-edge connection information.
 - Position information for each face (sufficient for unambiguous and continuous display; that is, shouldn't require optimization to get consistent connection positions)
 - Maybe: embedding info on another mesh?
 - Maybe: yarn ID info for knitting?


## Problems To Resolve

vertices represent corner-to-corner connections between faces that, in general, aren't represented by yarn.
(They also make it hard to translate `.knitout` directly to `.smobj`, since stacking stitches with xfer operations can lead to cases which don't make sense in the vertices-imply-connections case.)

Example:

```
>d>e>f>
 ^   ^
>a>b>c>
```
Stitches a,b,c sit in the course below stitches d,e,f; but stitch b is dropped before e was knit.
Unfortunately, the a^d and c^f connections imply that b's upper vertices correspond with e's lower vertices.

Face types often require yarn-reflected versions, which is an additional modeling burden; this could be made explicit in the format, at the price of additional loader complexity.
Knowledge of which faces are yarn-reflected versions is very useful when editing.

Faces might be parameterized by (e.g.) loop size, but this also makes corner-to-corner connection over-constrained.

Faces might also have different desired shapes depending on their neighbors, which implies a fair bit of modeling (or simulation?) burden.

(Low priority) As a text format, smobj may be slower to load/save than it could be, and may lose precision.
If there is a fast way to do better, it would be worth pursuing.

We want to make it easy to add extra data layers to the smobj (e.g., for correlation with photographs).
We could store these in separate files, or extra chunks, but certainly the text format helps with the extra-layers problem.

Proposed chunked format (variable-length faces are awkward!), edge-name format:
```
#basic geometry:
vXYZ NNNN
  (NNNN/12 vertices as little-endian f32 xyz tuples)
fSiz NNNN
  (NNNN face sizes, as u8)
fIdx NNNN
  (NNNN/4 indices into vertices as little-endian u32)
#edge type table:
etSt NNNN
  (NNNN bytes of utf8 edge type string data for edge type names)
etNm NNNN
  (NNNN/8 bytes of begin/end index pairs as little-endian u32 for edge type names)
#face type table:
ftSt NNNN
  (NNNN bytes of utf8 string data for face types)
ftNm NNNN
  (NNNN/8 bytes of begin/end index pairs as little-endian u32 for face type names)
  #need convention for edge directions; separate in/out types? seems odd to store per-exemplar.
ftSz NNNN
  (NNNN face type sizes, as u8)
ftEt NNNN
  (NNNN edge types per face, as u8 indicies into edge type table)
ftEd NNNN
  (NNNN edge directions per face, as i8 +1 / -1 / 0 for out, in, directionless)
#face types:
fTyp NNNN
  (NNNN face types, as u8 indices into face types table)
#might handle positions in images this way:
fUV0 NNNN
  (NNNN/8 uv coords per face as little-endian f32; NaN for no coordinates)
  #though it might be better to go full 3D coordinates.
```


## Some alternative ideas

Represent stitches by their frames and connection points; rendering done by warping based on connection points (instead of whole edges).
More awkward to display (or texture map) for sure.
Might generally be better for layout(?).

Keep faces but cut corners, giving each [e.g.] quad face four additional empty edges.

Keep faces but represent edge-edge connections explicitly; this means that code could at some point substitute in the cut-corners version.

Empty edges, in general, would complete the data structure nicely, in terms of consistency.
Do empty edges need special notation? (Probably yes -- empty edges either have no direction or perhaps can't be connected at all? Though allowing connections might be useful for visualization.)
