# smobj: a format for augmented-stitch-mesh-like things

Augmented Stitch Meshes represent knit structures and their dependencies by embedding yarns in the faces of a (tri/quad/pentagon/etc) mesh, and associating types and directions with each edge.

The `.smobj` format stores an augmented stitch mesh as follows:

```
#text format, like .obj with...
#vertices as X Y Z:
v 1.0 2.2 1.0
v 1.0 1.0 1.0
v 3.0 2.2 1.0
v 3.0 1.0 1.0
#faces as 1-based vertex index lists:
f 2 4 3 1
#and either:
#(1) explicit face names: (seems like a better format and less redundant)
T knit+
# ----
#or (2) one of 'tk', 'te', 'ts' to label edges where 0 = loopwise and 1= yarnwise
tk 1 0 1 0 
# and st to mark short-row faces, unused mostly
st 14
# 's' specifies starting face(s)
s f#
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
