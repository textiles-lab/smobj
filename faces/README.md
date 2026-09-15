# Stitch Mesh Face Libraries

This folder contains several face libraries for use in various stitch meshes.

## `illustration.sf` -- faces for knitting illustration

Faces included in the library:

### Increases (`make-*`)

Increases are point-up pentagonal faces with a loop entering from the bottom, yarns entering/exiting to the sides, and two loops exiting from the top.
The name indicates which of the existing loops is the "new" loop and how it is made, along with the direction of yarn travel.

Increases are named `make-`(`left`|`right`)`-`(`tuck`|`tuck-twist`|`split`)(`+`|`-`).
 - `left`|`right` indicates which loop is new (i.e., "not just knit through [or exactly] the old loop")
   - `left`: the left loop is the new one
   - `right`: the right loop is the new one
 - `tuck`|... indicates the method of making the new loop.
   - `tuck`: the new loop is tucked in the same direction as the yarn is already tavelling
   - `twisted-tuck`: the new loop is tucked on the front bed opposite the direction the yarn is already travelling. (I.e., a twisted tuck, with the direction of the twist such that the exiting yarn is behind the entering yarn.)
   - `split`: the new loop is split through the old loop (the other loop is just the old loop)
 - `+`|`-` indicates direction of yarn travel
   - `+` yarn enters from the left, exits to the right
   - `-` yarn enters from the right, exits to the left

*NOTE:* it is slightly odd that the `split` increase doesn't knit through the old loop after splitting through it, while the other increases all do this.

### Decreases (`decrease-*`)

Decreases are point-down pentogonal faces with two loops entering from the bottom, yarns entering/exiting from the sides, and a single loop (knit through both entering loops) exiting from the top.
The name indicates the stacking order of the loops and the direction of yarn travel.

Decreases are named `decrease-`(`left`|`right`)(`+`|`-`)
 - `left`|`right` indicates the lean of the decrease (i.e., which entering loop is stacked on top)
   - `left`: the right loop is stacked atop the left loop. (a "left-leaning" decrease)
   - `right`: the right loop is stacked atop the right loop. (a "right-leaning" decrease)
 - `+`|`-` direction of yarn travel
   - `+` yarn enters from the left, exits to the right
   - `-` yarn enters from the right, exits to the left

*NOTE:* in previous versions of this face library, the `-` versions of the decreases had the incorrect names (left vs right), such that changing from the `+` to the `-` version of a face would also change the loop stacking order.

### Turns (`turn-*`)

Faces for yarns turning at the end of a short row.
Diamond-shaped with yarns entering/exiting on two adjacent sides and loops entering/exiting on the other two sides.

Turns are named `turn`(`-tuck`)?(`)((`|`))(`)
 - `-tuck` indicates that the entering yarn tucks (in its direction of travel) a loop behind the entering loop. If missing, the yarn just turns around.
 - `)((`|`))(` indicates the shape
   - `)((`: the yarn enters via the lower left edge, exits via the upper left edge; the loop enters via the lower right edge, exits via the upper right edge.
   - `))(`: the yarn enters via the lower right edge, exits via the upper right face; the loop enters via the lower left edge, exits via the upper left edge.


### Edges (`edge*`)

Edges are pentagonal faces for connecting yarns between rows of knitting.

Edges are named `edge`(`(`|`)`).
 - `edge(` connects yarn from the lower right to the upper right face.
 - `edge)` connects yarn from the lower left to the upper left face.

### Starts (`start*`)

Starts are the first stitch made by an entering yarn. They are left-or-right-pointing triangular faces with loops entering from the bottom and exiting from the top, and yarns exiting to the left or right.

Starts are named `start`(`<`|`>`).
 - `start<` is a leftward-pointing triangle with yarn exiting to the right
 - `start>` is a rightward-pointing triangle with yarn exiting to the left

*NOTE:* the `.sf` format provides no way to start a yarn in the middle of a face, so the "fresh yarn" starts at the corner of the triangle.

### Ends (`end*`)

Ends are the last stitch made by an exiting yarn. They are left-or-right-pointing triangular faces with loops entering from the bottom and exiting from the top, and yarns entering to the left or right.

Ends are named `end`(`<`|`>`).
 - `end<` is a leftward-pointing triangle with yarn entering from the right
 - `end>` is a rightward-pointing triangle with yarn entering from the left

*NOTE:* the `.sf` format provides no way to end a yarn in the middle of a face, so the "exiting yarn" ends at the corner of the triangle.

### Ins / Outs (`in*`, `out*`)

Faces that just take a yarn in or out. Triangular with only one yarn in/out edge. The `+` variety is for yarns travelling right, the `-` variety is for yarns travelling left.

### Knits (`knit*`,`purl*`)

Rectangular stitches that front- or back-bed knit a stitch. Loops enter from the bottom and exit to the top; yarns enter from the left/right and exit to the right/left.

Knits are named (`knit`|`purl`)(`-`|`+`)
 - `knit`|`purl` indicates the direction the new loop is pulled through the old loop.
   - `knit`: the new loop is pulled forward (i.e., a front-bed knit)
   - `purl`: the new loop is pulled backward (i.e., a back-bed knit)
 - `+`|`-` direction of yarn travel
   - `+` yarn enters from the left, exits to the right
   - `-` yarn enters from the right, exits to the left

### Drops (`drop`)
Rectangular faces that drop a loop. A loop enters from the bottom, all other edges remain uncrossed.

*NOTE:* in the past, there were both `drop+` and `drop-` variations, but these don't actually differ unless you keep track of the directions of empty edges so were removed.

### Twisted Tucks (`tuck-twist*`)

Rectangular face holding an against-yarn-direction front-bed tuck (that is, a twisted tuck where the entering yarn passes in front of the exiting yarn).
The bottom edge is blank, the left and right edges connect to yarns, a loop exits from the top edge.

Twisted tucks are named `tuck-twist`(`+`|`-`)
 - `+`|`-` direction of yarn travel (also sets twist direction of tuck)
   - `+` yarn enters from the left, exits to the right
   - `-` yarn enters from the right, exits to the left

### Et Cetera

There are other faces in `illustration.sf` from various projects; in many cases there aren't full directional sets or they aren't quite modelled properly. As such, they are not documented here. (Eventually, though, they should be!)


## `knitout.sf` -- faces used by `knitout-to-smobj`

Basic faces for machine operations along with faces for yarn movement between operations.
Really not meant to be used other than as part of the `knitout-to-smobj` + `smobj-to-yarns` pipeline.

## `weave.sf` -- faces for weaving

Contains two square faces, `over` and `under`, modelling a weft crossing either over or under a warp.
