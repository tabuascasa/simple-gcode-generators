#!/usr/bin/python
# 18 Feb 2019 
# Minor edits

# based on the following:
# engrave-eds-meter-notice.py
# 16 Feb 2019
# I need different font scaling for some lines
# So repeat the initialisation of #1004 and #1005 with Xscale and Yscale
# ( Insert scale(xscale,yscale) before 'doText' calls to change scale)
# Added title to engrave constructor

# 15 Feb 2019
# engrave-lines3.1.py
# Based on engrave-lines Rev v2 21.06.2012 ArcEye
# Now works with python3
# JWB Jan 2019
# 
# Engrave a given string at X,Y in a chosen font and scaling
# See class 'engrave' below.


"""
    engraver.py: Copyright (C) 2019 Joe Brown, is based on work by several
    other authors. Previously copyrighted material is mentioned below.

    engraver.py will produce CNC gcode to engrave multiple text lines.
    The module is fully python3 compatible.
    The module is presented as classes, the principal of which is 'engrave'.
    An example of use is presented in 'main()' at the bottom of this file.

    The original GUI code has been removed, as the output code can be tested
    properly on simulations using linuxcnc et al, very quickly.
    The original cmd line code has been removed as being inappropriate for 
    including the module into current work.
    Where possible, original variable names have been retained.
    Most of the features of the original have been retained, but not all
    of these have been re-tested by me.

    Features in this version:
     * Gcode can be output to either stdout or a file, or both.
     * Individual fonts can be used on each separate text line.
     * Individual lines can be scaled in both X and Y.
     * X and Y start positions for individual line text can be specified,
       including a text centreing option.
     * Output files contain both a given title and date-time stamp.
     * Pre and Postamble code strings can be easily edited to suit
       individual CNC software and/or machines

    Original Copyright/Use Notice:
    Based upon code from engrave-11.py
    Copyright (C) <2008>  <Lawrence Glaister> <ve7it at shaw dot ca>
    based on work by John Thornton  -- GUI framwork from arcbuddy.py
    Ben Lipkowitz  (fenn)-- cxf2cnc.py v0.5 font parsing code

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    Rev v2 21.06.2012 ArcEye
"""

from math import *
import re
import sys
from datetime import datetime

#=======================================================================
class Character:

    def __init__(self, key):
        self.key = key
        self.stroke_list = []

    def __repr__(self):
        return "%s" % (self.stroke_list)

    def get_xmax(self):
        try: 
            return max([s.xmax for s in self.stroke_list[:]])
        except ValueError: 
            return 0

    def get_ymax(self):
        try: 
            return max([s.ymax for s in self.stroke_list[:]])
        except ValueError: 
            return 0

#=======================================================================
class Line:

    def __init__(self, coords):
        self.xstart, self.ystart, self.xend, self.yend = coords
        self.ymax = max(self.ystart, self.yend)
        self.xmax = max(self.xstart, self.xend)

    def __repr__(self):
        return "Line([%s, %s, %s, %s])" % (self.xstart, self.ystart, self.xend, self.yend)

#=======================================================================
class engrave:

    SafeZ = 5.0
    Depth = 0.1
    XScale =  1
    YScale =  1
    CSpaceP = 25
    WSpaceP =  100
    Angle = 0
    Mirror = False
    Flip = False
    Postamble = ["G40", "( Home the tool )", "G53 G01 X0 Y0 Z{:03.3f} M5".format(SafeZ), "M30"]
    outtostd = True
    outtofile = False

    visit = 0
    gcode = []

    # change this if you want to default another font
    default_fontfile = "/Users/joebrown/Documents/cxf-fonts/romanc.cxf"
    current_fontfile = ""
    font = None
    font_line_height = 0.0
    font_word_space =  0.0
    font_char_space = 0.0

    def __init__(self, title, jobposx=0, jobposy=0, safeZ=5, depth=0.1, xscale=1, yscale=1, cspacep=25, wspacep=100, angle=0, mirror=False, flip=False):
        self.SafeZ = safeZ
        self.Depth = depth
        self.XScale =  xscale
        self.YScale =  yscale
        self.CSpaceP = cspacep
        self.WSpaceP =  wspacep
        self.Angle = angle
        self.Mirror = mirror
        self.Flip = flip
        self.jobx = jobposx
        self.joby = jobposy
        self.setjobpos(jobposx,jobposy)
        self.fws_scale = 0.5
        self.title = title
        
        datestring = datetime.strftime(datetime.now(), '%Y-%m-%d-%H-%M-%S')
        self.filestring = "engrave_" + title + "_" + datestring + ".ngc"
        self.linenum = 1

    def outopts(self, tostd=True, tofile=False):
        self.outtostd = tostd
        self.outtofile = tofile

    def doPostamble(self, gcode):
        for amble in self.Postamble:
            gcode.append(amble)

    def setjobpos(self, xoff, yoff):
        self.jobpos = "X{} Y{}".format("{:03.3f}".format(self.jobx), "{:03.3f}".format(self.joby))
        self.Preamble = "G40 G10 L2 P1 " + self.jobpos + " G54 G21 G90 G94 F50"

    def closefile(self):
        self.outfile.close()

    def depth(self, cutdepth):
        self.Depth = cutdepth

    def safeZ(self, safez):
        self.SafeZ = safez

    def angle(self, newangle):
        self.Angle = newangle

    def scale(self, xscale, yscale):
        self.XScale = xscale
        self.YScale = yscale

    def flip(self, flip):
        self.Flip = flip

    def mirror(self, mirror):
        self.Mirror = mirror

    def spaces(self, cspace, wspace):
        self.CSpaceP = cspace
        self.WSpaceP =  wspace        

    def preamble(self,preambstr):
        self.Preamble = preambstr

    def postamble(self,postambstr):
        self.Postamble = postambstr

    def defaultFont(self, fontpath):
        self.default_font = fontpath

    def sanitize(self, string):
        retval = ''
        good=' ~!@#$%^&*_+=-{}[]|\:;"<>,./?'
        for char in string:
            if char.isalnum() or good.find(char) != -1:
                retval += char
            else: 
                retval += ( ' 0x%02X ' %ord(char))
        return retval

    #=======================================================================
    # This routine parses the .cxf font file and builds a font dictionary of
    # line segment strokes required to cut each character.
    # Arcs (only used in some fonts) are converted to a number of line
    # segemnts based on the angular length of the arc. Since the idea of
    # this font description is to make it support independant x and y scaling,
    # we can not use native arcs in the gcode.
    #=======================================================================
    def parse(self, ffile):
        font = {}
        key = None
        num_cmds = 0
        line_num = 0
        lines = ffile.readlines()
        for text in lines:
            #format for a typical letter (lowercase r):
            ##comment, with a blank line after it
            #
            #[r] 3
            #L 0,0,0,6
            #L 0,6,2,6
            #A 2,5,1,0,90
            #
            line_num += 1
            end_char = re.match('^$', text) #blank line
            if end_char and key: #save the character to our dictionary
                font[key] = Character(key)
                font[key].stroke_list = stroke_list
                font[key].xmax = xmax
                if (num_cmds != cmds_read):
                    print ("(warning: discrepancy in number of commands %s, line %s, %s != %s )" % (fontfile, line_num, num_cmds, cmds_read))

            new_cmd = re.match('^\[(.*)\]\s(\d+)', text)
            if new_cmd: #new character
                key = new_cmd.group(1)
                num_cmds = int(new_cmd.group(2)) #for debug
                cmds_read = 0
                stroke_list = []
                xmax, ymax = 0, 0

            line_cmd = re.match('^L (.*)', text)
            if line_cmd:
                cmds_read += 1
                coords = line_cmd.group(1)
                coords = [float(n) for n in coords.split(',')]
                stroke_list += [Line(coords)]
                xmax = max(xmax, coords[0], coords[2])

            arc_cmd = re.match('^A (.*)', text)
            if arc_cmd:
                cmds_read += 1
                coords = arc_cmd.group(1)
                coords = [float(n) for n in coords.split(',')]
                xcenter, ycenter, radius, start_angle, end_angle = coords
                # since font defn has arcs as ccw, we need some font foo
                if ( end_angle < start_angle ):
                    start_angle -= 360.0
                # approximate arc with line seg every 20 degrees
                segs = int((end_angle - start_angle) / 20) + 1
                angleincr = (end_angle - start_angle)/segs
                xstart = cos(start_angle * pi/180) * radius + xcenter
                ystart = sin(start_angle * pi/180) * radius + ycenter
                angle = start_angle
                for i in range(segs):
                    angle += angleincr
                    xend = cos(angle * pi/180) * radius + xcenter
                    yend = sin(angle * pi/180) * radius + ycenter
                    coords = [xstart,ystart,xend,yend]
                    stroke_list += [Line(coords)]
                    xmax = max(xmax, coords[0], coords[2])
                    ymax = max(ymax, coords[1], coords[3])
                    xstart = xend
                    ystart = yend

        self.font_line_height = max(font[key].get_ymax() for key in font)
        self.font_word_space =  max(font[key].get_xmax() for key in font) * (self.WSpaceP / 100.0)
        self.font_word_space = self.font_word_space * self.fws_scale

        self.font_char_space = self.font_word_space * (self.CSpaceP / 100.0)

        return font

    def doText(self, strtoengrave, x, y, centreX=None, ff=None):
        self.gcode = [] # clear last result
        lenx = 0 # length of string in mm
        givenX = x

        # setup font file
        changed = False
        if ff == "current": # no change requested
            if self.current_fontfile == "":
                self.current_fontfile = self.default_fontfile
                changed = True

        elif ff == None: # use default font
            if self.current_fontfile != self.default_fontfile:
                self.current_fontfile = self.default_fontfile
                changed = True
        
        elif ff != self.current_fontfile: # font not already set up, or different
            self.current_fontfile = ff
            changed = True

        # NOTE: use same font if NOT changed        
        if changed:
            # NOTE ENCODING!! latin-1 allows ALL 256 characters possible in a byte - stops annoying UnicodeDecode errors
            ffile = open(self.current_fontfile, encoding="latin-1")
            self.font = self.parse(ffile)          # build stroke lists from font file
            ffile.close()

        # process given x coordinate if centreX is true
        # JWB Note derivative x is an approximation only - my knowledge of fonts is sketchy at best
        if centreX == True:
            # x given is for the centre of the word, not the beginning
            for char in strtoengrave:
                if char == ' ':
                    lenx += self.font_word_space
                else:
                    char_width = self.font[char].get_xmax() * self.XScale
                    lenx +=  char_width  # + self.font_char_space
            lstr = len(strtoengrave)
            if lstr % 2 == 0:
                avg = lenx / lstr
                lenx += avg 
            x = x - lenx/2 # adjust start of word to x - half length of word in mm
                
            if x < 0:
                print("Centreing " + strtoengrave + " results in NEGATIVE x coordinate, exiting")
                sys.exit()

        engravingblah = []
        engravingblah.append("( Engraving: '{}' at X{:.3f}, Y{:.3f})".format(strtoengrave, x, y))
        if centreX == True:
           engravingblah.append("( '{}' is centred on X{:.3f} )".format(strtoengrave, givenX))
        engravingblah.append("( using font file: {} )".format(self.current_fontfile))

        if self.visit != 0: # subsequent call to routine
            self.gcode.append("(===================================================================)")
            self.gcode.extend(engravingblah)
            self.gcode.append("#1002 = {:.4f}  ( X Start )".format(x))
            self.gcode.append("#1003 = {:.4f}  ( Y Start )".format(y))
            self.gcode.append("#1004 = {:.4f}  ( X Scale )".format(self.XScale))
            self.gcode.append("#1005 = {:.4f}  ( Y Scale )".format(self.YScale))
            self.gcode.append("(===================================================================)")
            
        else: # 1st call to routine since class instantiated
            if self.outtofile:
                self.outfile = open(self.filestring, 'w')
                self.gcode.append("( File: {} )".format(self.filestring))

            self.gcode.append('( Code generated using engraver.py )')
            self.gcode.append('( Authors: ArcEye, Lawrence Glaister, and Joe Brown)')
            self.gcode.append('( This version: Joe Brown Jan 2019. )')
            self.gcode.extend(engravingblah)
            # write out subroutine for rotation logic just once at head of gcode
            self.gcode.append("(===================================================================)")
            self.gcode.append("(Subroutine to handle x,y rotation about 0,0)")
            self.gcode.append("(input x,y get scaled, rotated then offset )")
            self.gcode.append("( [#1 = 0 or 1 for a G0 or G1 type of move], [#2=x], [#3=y])")
            self.gcode.append("o9000 sub")
            self.gcode.append("  #28 = [#2 * #1004]  ( scaled x )")
            self.gcode.append("  #29 = [#3 * #1005]  ( scaled y )")
            self.gcode.append("  #30 = [SQRT[#28 * #28 + #29 * #29 ]]   ( dist from 0 to x,y )")
            self.gcode.append("  #31 = [ATAN[#29]/[#28]]                ( direction to  x,y )")
            self.gcode.append("  #32 = [#30 * cos[#31 + #1006]]     ( rotated x )")
            self.gcode.append("  #33 = [#30 * sin[#31 + #1006]]     ( rotated y )")
            self.gcode.append("  o9010 if [#1 LT 0.5]" )
            self.gcode.append("    G00 X[#32+#1002] Y[#33+#1003]")
            self.gcode.append("  o9010 else")
            self.gcode.append("    G01 X[#32+#1002] Y[#33+#1003]")
            self.gcode.append("  o9010 endif")
            self.gcode.append("o9000 endsub")
            self.gcode.append("(===================================================================)")
        
            self.gcode.append("#1000 = {:.4f}".format(self.SafeZ))
            self.gcode.append("#1001 = {:.4f}  ( Engraving Depth Z )".format(self.Depth))
            self.gcode.append("#1002 = {:.4f}  ( X Start )".format(x))
            self.gcode.append("#1003 = {:.4f}  ( Y Start )".format(y))
            self.gcode.append("#1004 = {:.4f}  ( X Scale )".format(self.XScale))
            self.gcode.append("#1005 = {:.4f}  ( Y Scale )".format(self.YScale))
            self.gcode.append("#1006 = {:.4f}  ( Angle )".format(self.Angle))
            self.gcode.append(self.Preamble)
            
        self.gcode.append( 'G0 Z#1000')

        xoffset = 0 # distance along raw string in font units
        word = ""
        oldx = oldy = -99990.0      

        for char in strtoengrave:
            if char == ' ':
                xoffset += self.font_word_space
                self.gcode.append("( Completed word: '{}' )".format(word))
                word = ""
                continue
                          
            try:
                word += char
                self.gcode.append("(character '%s')" % self.sanitize(char))

                first_stroke = True
                for stroke in self.font[char].stroke_list:
                    dx = oldx - stroke.xstart
                    dy = oldy - stroke.ystart
                    dist = sqrt(dx * dx + dy * dy)

                    x1 = stroke.xstart + xoffset
                    y1 = stroke.ystart
                    if self.Mirror:
                        x1 = -x1
                    if self.Flip:
                        y1 = -y1

                    # check and see if we need to move to a new discontinuous start point
                    if (dist > 0.001) or first_stroke:
                        first_stroke = False
                        #lift engraver, rapid to start of stroke, drop tool
                        self.gcode.append("G0 Z#1000")
                        self.gcode.append("o9000 call [0] [{:.4f}] [{:.4f}]".format(x1,y1))
                        self.gcode.append("G1 Z#1001")

                    x2 = stroke.xend + xoffset
                    y2 = stroke.yend
                    if self.Mirror:
                        x2 = -x2
                    if self.Flip:
                        y2 = -y2
                    self.gcode.append("o9000 call [1] [{:.4f}] [{:.4f}]".format(x2,y2))
                    oldx, oldy = stroke.xend, stroke.yend

                # move over for next character
                char_width = self.font[char].get_xmax()
                xoffset += self.font_char_space + char_width

            except KeyError:
               self.gcode.append("(warning: character '0x%02X' not found in font defn)" % ord(char))

            self.gcode.append("")       # blank line after every char block

        self.gcode.append("( Completed word: '{}' )".format(word))
        self.gcode.append( 'G0 Z#1000')     # final engraver up
        # bump instance calls
        self.visit += 1
        self.outputGCode()

    def outputGCode(self):
        if self.outtofile:
            for line in self.gcode:
                # include line numbers in file
                linenumstr = "N{}".format(self.linenum)
                self.outfile.write(linenumstr + " " + line +'\n')
                self.linenum += 1

        if self.outtostd:      
            for line in self.gcode:
                sys.stdout.write(line+'\n')

    def finish(self):
        self.gcode = []
        self.doPostamble(self.gcode)
        self.outputGCode()
        if self.outtofile:
            self.closefile()

#===============================================================================================================

def main():

    Xcentre = 150.00   # centre the text line(s) at 150mm
    #                  title - will be part of output filename
    gencode = engrave("EDP_Meter_Notice")
    gencode.outopts(False,True) # no output to stdout, output to file
    gencode.depth(-0.1) # engraving depth
    
    #====================================================================================================
    # simple test - produce gcode for text engraving small panel indicating location of electricity meter
    #
    #====================================================================================================
    gencode.scale(3.0,3.0) # 3 X the font size
    #               string X,       Y,     X is centre of string, font
    gencode.doText("EDP",  Xcentre, 70.00, True) # no font: default font (romanc)
    gencode.scale(2.0,2.0) # 2 X the font size
    gencode.doText("Medidor de eletricidade", Xcentre,  45, True, "current") # current: current font
    gencode.doText("esta por tras deste obturador", Xcentre,  20.00, True, "current")               
    gencode.finish()
            
#===============================================================================================
            
if __name__ == "__main__":
	    main()
