#path to cp2k executable
PATH_CP2K=$1

# The xml file stores all the documentation data (i.e. structure of the code etc...)
XML_FILE="cp2k_input.xml"

# The XSL file is a stylesheet file that tells your browser how to present the HTML code.
XSL_FILE="cp2k_input.xsl"

# The java code for the saxon executable
SAXON_EXE="saxon-he-10.3.jar"

# The directory which we save everything to
OUT_DIR="Docs"




# Just do some checks 
if [ "$PATH_CP2K" == "" ]
then
	echo "Please give the cp2k path and the executable: e.g. path/cp2k.sopt."
        exit
fi

if [ "$SAXON_EXE" == "" ]
then
	echo "I can't find the variable 'SAXON_EXE'."
	echo "Please declare it at the top of this script e.g.:"
	echo "'SAXON_EXE="saxon-he-10.3.jar""
fi

if ! [ -f "$SAXON_EXE" ]
then
	echo "I can't find the saxon executable. Please check for typos!"
	echo "SAXON_EXE = '$SAXON_EXE'"
fi

if [ "$XML_FILE" == "" ]
then
	echo "I can't find the XML file variable -please link it in the code by declaring a variable at the top named 'XML_FILE' e.g. 'XML_FILE="cp2k_input.xml"'"
	exit
fi

if [ "$XSL_FILE" == "" ]
then
	echo "I can't find the XSL file variable -please link it in the code by declaring a variable at the top named 'XSL_FILE' e.g. 'XSL_FILE="cp2k_input.xml"'"
	exit
fi

if ! [ -f "$XSL_FILE" ]
then
	echo "The specified XSL file doesn't exist. Please check for typos!"
	echo "XSL_FILE = $XSL_FILE"
	exit
fi

#create xml file
echo "path" $PATH_CP2K
"$PATH_CP2K" --xml 
echo "xml file created"


if ! [ -f "$XML_FILE" ]
then
	echo "The specified XML file doesn't exist. Please check for typos!"
	echo "XML_FILE = $XML_FILE"
	exit
fi

# Create our directory and fill it
if [ "$OUT_DIR" == "" ]
then
	OUT_DIR="Docs"
fi
if [ -d "$OUT_DIR" ]
then
	rm -rf $OUT_DIR
fi
mkdir $OUT_DIR
cp $SAXON_EXE $XML_FILE $XSL_FILE $OUT_DIR
cd $OUT_DIR


# Now run the saxon code
java -jar $SAXON_EXE -s:$XML_FILE -xsl:$XSL_FILE -o:CP2K_Documentation.html
rm $SAXON_EXE $XML_FILE $XSL_FILE

