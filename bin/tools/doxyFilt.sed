# -----------------------------------------------------------------------------
# Script
#     doxyFilt.sed
#
# Description
#     Transform human-readable tags such as 'Description' into the Doxygen
#     equivalent
# -----------------------------------------------------------------------------

/^License/,/\*\//{
/^License/,/MA 0211.-130. USA/{
s?^License.*?\*\/\
\/\*! @file %filePath%\
<b>Original source file</b> <a href="%filePath%">%fileName%</a>\
?
/^    /d
}

# remove entry
/^Primitive *$/{
N
N
d
}

# remove entry
/^Implementation *$/{
N
N
d
}

# remove entry
/^Application *$/{
N
N
d
}

# remove entry
/^Type *$/{
N
N
d
}

# remove entry
/^Global *$/{
N
N
d
}


# Class
#     Foam::className
# =>
# @class Foam::className
#
/^Class *$/,/^[^ ]/{
/^Class/d
s/^    /@class /
}


# Namespace
#     namespaceName
# =>
# @namespace namespaceName
#
/^Namespace *$/,/^[^ ]/{
/^Namespace/d
s/^    /@namespace /
}


# Typedef
#     Foam::def
# =>
# @class Foam::def
# This is not strictly correct, but makes it easier to find the typedefs
/^Typedef *$/,/^[^ ]/{
/^Typedef/d
s/^    /@class /
}


# add anchor and use @brief
# the first paragraph will be 'brief' and the others 'detail'
/^Description *$/,/^[^ ]/{
/^Description/c\
<a class="anchor" name="Description"></a>\
@brief
s/^    //
}

/^Usage *$/,/^[^ ]/{
/^Usage/c\
@par Usage
s/^    //
}


/^See *Also *$/,/^[^ ]/{
/^See *Also/c\
@see
s/^    //
}

/^Author *$/,/^[^ ]/{
/^Author/c\
@author
s/^    //
}

/^Note *$/,/^[^ ]/{
/^Note/c\
@note
s/^    //
}


/^To[Dd]o *$/,/^[^ ]/{
/^To[Dd]o/c\
@todo
s/^    //
}

/^Warning *$/,/^[^ ]/{
/^Warning/c\
@warning
s/^    //
}

/^Deprecated *$/,/^[^ ]/{
/^Deprecated/c\
@deprecated
s/^    //
}

/SourceFiles/,/^[ ]*$/{
s?SourceFiles?@par Source files\
<ul>\
  <li><a href="%filePath%">%fileName%</a></li>?
s?^[ ]*$?</ul>\
?
s? *\([a-zA-Z0-9]*\.[a-zA-Z]*\)?  <li><a href="%dirName%/\1">\1</a></li>?
}

s/.*\*\//\*\//

}

# -----------------------------------------------------------------------------
