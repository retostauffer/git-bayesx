
//package gov.sandia.postscript;

import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import java.awt.event.*;

/**
 * PSGr is a Graphics subclass that images to PostScript.
 * (C) 1996 E.J. Friedman-Hill and Sandia National Labs
 * @version 	2.1
 * @author 	Ernest Friedman-Hill
 * @author      ejfried@ca.sandia.gov
 * @author      http://herzberg.ca.sandia.gov
 */

public abstract class PSGrBase extends java.awt.Graphics
{
//  public static final boolean DEBUG = true;
  public static final boolean DEBUG = false;    
/*
  protected final static int PAGEHEIGHT = 792;
  protected final static int PAGEWIDTH = 612;
  protected final static int XOFFSET = 30;
  protected final static int YOFFSET = 30;
*/
  protected final static int PAGEHEIGHT = 842;
  protected final static int PAGEWIDTH = 596;
  protected final static int XOFFSET = 0;
  protected final static int YOFFSET = 0;
  
/**
     hexadecimal digits
  */
  
  protected final static char hd[] = {'0', '1', '2', '3', '4', '5', '6', '7',
                                      '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
  
  /**
     number of chars in a full row of pixel data
  */
  
  protected final static int charsPerRow = 12*6;

  
  /**
     Output stream where postscript goes
  */

  protected PrintWriter os;

  /**
     The current color
  */

  protected Color clr = Color.black;


  /**
     The background color of the current widget.
     It's up to the client software to set this correctly!
  */

  protected Color backClr = Color.white;

  /**
     The current font
  */

  protected Font font = new Font("Helvetica",Font.PLAIN,12);

  protected Rectangle clippingRect = new Rectangle(0,0,PAGEWIDTH,PAGEHEIGHT);

  /**
   * Constructs a new PSGr Object that write to stdout. Unlike regular Graphics objects,
   * PSGr contexts can be created directly.
   */

  public PSGrBase()
  {
    this(new OutputStreamWriter(System.out), true);
  }

  /**
   * Constructs a new PSGr Object. Unlike regular Graphics objects,
   * PSGr contexts can be created directly.
   * @param o Output stream for PostScript output
   */

  public PSGrBase(Writer w)
  {
    this(w, true);
  }

  /**
   * Constructs a new PSGr Object. Unlike regular Graphics objects,
   * PSGr contexts can be created directly.
   * @param o Output stream for PostScript output
   */

  public PSGrBase(Writer w, boolean emitProlog)
  {
    setOutput(w);
    
    if (emitProlog)
      emitProlog();
  }

  /**
   * Creates a new PSGr Object that is a copy of the original PSGr Object.
   * Not implemented; throws RuntimeException
   */
  public Graphics create()
  {
    throw new RuntimeException("create() not implemented");
  }

  /**
   * Change the Writer this context's output goes to.
   */

  public void setOutput(Writer w)
  {
    if (! (w instanceof PrintWriter))
      os = new PrintWriter(w, true);
    else
      os = (PrintWriter) w;    
  }

  
  /**
   * Translates the specified parameters into the origin of
   * the graphics context. All subsequent
   * operations on this graphics context will be relative to this origin.
   * @param x the x coordinate
   * @param y the y coordinate
   * @see #scale
   */

  public void translate(int x, int y)
  {
    os.print(x);
    os.print(" ");
    os.print(y);
//    os.println(" translate");
    os.print(" T ");    
  }
  
  /**
   * Scales the graphics context. All subsequent operations on this
   * graphics context will be affected.
   * @param sx the scaled x coordinate
   * @param sy the scaled y coordinate
   * @see #translate
   */
  protected void scale(double x, double y)
  {
    os.print(x);
    os.print(" ");
    os.print(y);
//    os.println(" scale");
    os.println(" sc ");    
  }

  protected void lineto(int x, int y)
  {
    os.print(x);
    os.print(" ");
    os.print(y);
//    os.println(" lineto");
    os.print(" L ");    
  }

  protected void moveto(int x, int y)
  {
    os.print(x);
    os.print(" ");
    os.print(y);
//    os.println(" moveto");
    os.print(" M ");
  }

  protected void lineto(double x, double y)
  {
    os.print(x);
    os.print(" ");
    os.print(y);
//    os.println(" lineto");
    os.print(" L ");        
  }

  protected void moveto(double x, double y)
  {
    os.print(x);
    os.print(" ");
    os.print(y);
//    os.println(" moveto");
    os.print(" M ");
  }
  

  /**
   * Gets the current color.
   * @see #setColor
   */
  public Color getColor()
  {
    if (DEBUG) diagnostic("getColor()");
    return clr;
  }

  /**
   * Gets the current color.
   * @see #setColor
   */
  public void setBackground(Color c)
  {
    if (DEBUG) diagnostic("setBackground(" + c + ")");
    backClr = c;;
  }
  

  /**
   * Sets the current color to the specified color. All subsequent graphics operations
   * will use this specified color.
   * @param c the color to be set
   * @see Color
   * @see #getColor
   */

  public void setColor(Color c)
  {
    if (DEBUG) diagnostic("setColor(" + c + ")");
    if (c != null)
      clr = c;
    os.print(clr.getRed()/255.0);
    os.print(" ");
    os.print(clr.getGreen()/255.0);
    os.print(" ");
    os.print(clr.getBlue()/255.0);
    os.println(" setrgbcolor");
  }

  /**
   * Sets the default paint mode to overwrite the destination with the
   * current color. PostScript has only paint mode.
   */
  public void setPaintMode()
  {
    if (DEBUG) diagnostic("setPaintMode()");
  }


  /**
   * Sets the paint mode to alternate between the current color
   * and the new specified color. PostScript does not support XOR mode.
   * @param c1 the second color
   */
  public void setXORMode(Color c1)
  {
    if (DEBUG) diagnostic("setXORMode(" + c1 + ")");
  }

  /**
   * Gets the current font.
   * @see #setFont
   */
  public Font getFont()
  {
    if (DEBUG) diagnostic("getFont()");
    return font;
  }

  /**
   * Sets the font for all subsequent text-drawing operations.
   * @param font the specified font
   * @see Font
   * @see #getFont
   * @see #drawString
   * @see #drawBytes
   * @see #drawChars
   */
  public void setFont(Font f)
  {
    if (DEBUG) diagnostic("setFont(" + f + ")");

    if (f != null) {
      this.font = f;
      String fontName = font.getName();
      int fontStyle = font.getStyle();

      if (fontName.startsWith("Times"))
        {
          switch (fontStyle)
            {
            case Font.PLAIN:
              fontName = "Times-Roman"; break;
            case Font.BOLD:
              fontName = "Times-Bold"; break;
            case Font.ITALIC:
              fontName = "Times-Italic"; break;
            case (Font.ITALIC + Font.BOLD):
              fontName += "Times-BoldItalic"; break;
            }
        }     
      else     
        {
          switch (fontStyle)
            {
            case Font.PLAIN:
              break;
            case Font.BOLD:
              fontName += "-Bold"; break;
            case Font.ITALIC:
              fontName += "-Oblique"; break;
            case (Font.ITALIC + Font.BOLD):
              fontName += "-BoldOblique"; break;
            }
        }
      
      os.println("/" + fontName + " findfont");
      os.print(font.getSize());
      os.println(" scalefont setfont");
    }
  }

  /**
   * Gets the current font metrics.
   * @see #getFont
   */
  public FontMetrics getFontMetrics()
  {
    if (DEBUG) diagnostic("getFontMetrics()");
    return getFontMetrics(getFont());
  }

  /**
   * Gets the current font metrics for the specified font.
   * @param f the specified font
   * @see #getFont
   * @see #getFontMetrics
   */
  public FontMetrics getFontMetrics(Font f)
  {
    if (DEBUG) diagnostic("getFontMetrics(" + f + ")");
    return Toolkit.getDefaultToolkit().getFontMetrics(f);
  }


  /** 
   * Returns the bounding rectangle of the current clipping area.
   * @see #clipRect
   * @see #getClipBounds
   * @deprecated
   */
  public Rectangle getClipRect()
  {
    if (DEBUG) diagnostic("getClipRect()");
    return clippingRect;
  }

  /** 
   * Returns the bounding rectangle of the current clipping area.
   * @see #clipRect
   * @see #getClipBounds
   */
  public Shape getClip()
  {
    if (DEBUG) diagnostic("getClip()");
    return clippingRect;
  }

  /** 
   * Sets the clipping region using the shape's bounding rectangle.
   * @param s A shape to set the clipping rectangle to
   */
  public void setClip(Shape s)
  {
    if (DEBUG) diagnostic("setClip(" + s + ")");
    Rectangle r = s.getBounds();
    setClip(r.x, r.y, r.width, r.height);
  }


  /** 
   * Returns the bounding rectangle of the current clipping area.
   * @see #clipRect
   * @see #getClipRect
   * @return the bounding rectangle
   */
  public Rectangle getClipBounds()
  {
    if (DEBUG) diagnostic("getClipBounds()");
    return clippingRect;
  }

  /** 
   * Clips to a rectangle. The resulting clipping area is the
   * intersection of the current clipping area and the specified
   * rectangle. Graphic operations have no effect outside of the
   * clipping area.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @see #getClipRect
   */
  public void setClip(int x, int y, int width, int height)
  {
    if (DEBUG) diagnostic("setClip("+ x + ", " +y + ", " +width + ", " + height +")");
    y = transformY(y);
    clippingRect = new Rectangle(x,y,width,height);
    os.println("initclip");

    moveto(x,y);
    lineto(x + width, y);
    lineto(x + width, y - height);
    lineto(x, y - height);
    os.println("closepath eoclip newpath");
  }

  /** 
   * @deprecated
   */
  public void clipRect(int x, int y, int width, int height)
  {
    if (DEBUG) diagnostic("clipRect(" + x + ", " +y + ", " +width + ", " +height +")");
    setClip(x, y, width, height);
  }

  /**
   * Copies an area of the screen.
   * @param x the x-coordinate of the source
   * @param y the y-coordinate of the source
   * @param width the width
   * @param height the height
   * @param dx the horizontal distance
   * @param dy the vertical distance
   * Note: copyArea not supported by PostScript
   */
  public void copyArea(int x, int y, int width, int height, int dx, int dy)
  {
    if (DEBUG) diagnostic("copyArea(" + x + ", " +y + ", " +width + ", " +height + ", " +dx + ", " +dy +")");
    throw new RuntimeException("copyArea not supported");
  }

  /** 
   * Draws a line between the coordinates (x1,y1) and (x2,y2). The line is drawn
   * below and to the left of the logical coordinates.
   * @param x1 the first point's x coordinate
   * @param y1 the first point's y coordinate
   * @param x2 the second point's x coordinate
   * @param y2 the second point's y coordinate
   */
  public void drawLine(int x1, int y1, int x2, int y2) 
  {
    if (DEBUG) diagnostic("drawLine(" + x1 + ", " +y1 + ", " +x2 + ", " +y2 +")");
    y1 = transformY(y1);
    y2 = transformY(y2);
    moveto(x1, y1);
    lineto(x2, y2);
    stroke(false);
  }

  public void drawLine(double x1, double y1, double x2, double y2) 
  {
    if (DEBUG) diagnostic("drawLine(" + x1 + ", " +y1 + ", " +x2 + ", " +y2 +")");
    y1 = transformY(y1);
    y2 = transformY(y2);
    moveto(x1, y1);
    lineto(x2, y2);
    stroke(false);
  }  
  
  protected void doRect(int x, int y, int width, int height, boolean fill) 
  {
    if (DEBUG) diagnostic("doRect(" + x + ", " +y + ", " +width + ", " +height + ", " +fill+")");
    y = transformY(y);
    moveto(x, y);
    lineto(x + width, y);
    lineto(x + width, y - height);
    lineto(x, y-height);
    lineto(x, y);
    stroke(fill);      
  }

  /** 
   * Fills the specified rectangle with the current color. 
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @see #drawRect
   * @see #clearRect
   */
  public void fillRect(int x, int y, int width, int height) 
  {
    if (DEBUG) diagnostic("fillRect(" + x + ", " +y + ", " +width + ", " +height+")");
    doRect(x,y,width,height,true);
  }

  /** 
   * Draws the outline of the specified rectangle using the current color.
   * Use drawRect(x, y, width-1, height-1) to draw the outline inside the specified
   * rectangle.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @see #fillRect
   * @see #clearRect
   */
  public void drawRect(int x, int y, int width, int height) 
  {   
    if (DEBUG) diagnostic("drawRect(" + x + ", " +y + ", " +width + ", " +height+")");
    doRect(x,y,width,height,false);
  }
  
  /** 
   * Clears the specified rectangle by filling it with the current background color
   * of the current drawing surface.
   * Which drawing surface it selects depends on how the graphics context
   * was created.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @see #fillRect
   * @see #drawRect
   */
  public void clearRect(int x, int y, int width, int height) 
  {
    if (DEBUG) diagnostic("clearRect(" + x + ", " +y + ", " +width + ", " +height+")");
    gsave();
    Color c = getColor();
    setColor(backClr);
    doRect(x,y,width,height, true);
    setColor(c);
    grestore();
  }


  protected void doRoundRect(int x, int y, int width, int height,
                           int arcWidth, int arcHeight, boolean fill) 
  {
    if (DEBUG) diagnostic("doRoundRect(" + x + ", " +y + ", " +width + ", " +height
                          +arcWidth + ", " +arcHeight + ", " + fill +")");

    y = transformY(y);

    gsave();

    // This value is OK if the two arc dimensions are the same
    int arcDim = arcHeight/2;

    translate(x, y);
        
    if (arcHeight != arcWidth)
      {

        // Postscript can't draw elliptical arcs directly. We have to scale a
        // rectangle with circular arcs to get what we want.

        if (arcHeight > arcWidth)
          {
            double ratio = (double) arcHeight / (double) arcWidth;
            scale(1.0, ratio);
            height = (int) ((double) height / ratio);
            arcDim = arcWidth/2;
          }
        else
          {
            double ratio = (double) arcWidth / (double) arcHeight;
            scale(ratio, 1.0);
            width = (int) ((double) width / ratio);
            arcDim = arcHeight/2;
          }

      }

    // Draw the actual rectangle
    os.println("0 setlinewidth");
    moveto(arcDim, 0);    
    arcTo(width, 0, width, -height, arcDim);
    arcTo(width, -height, 0, -height, arcDim);
    arcTo(0,-height, 0, 0, arcDim);
    arcTo(0, 0, width, 0, arcDim);

    stroke(fill);
    os.println("1 setlinewidth");

    grestore();
  }

  protected void stroke(boolean fill)
  {
    if (fill) 
      {
        gsave();
//        os.println("eofill");
        os.print(" E ");
        grestore();
      }
    
//    os.println("stroke");
    os.println(" S ");    
  }

  protected void arcTo(int x1, int y1, int x2, int y2, int dim)
  {
    os.print(x1);
    os.print(" ");
    os.print(y1);
    os.print(" ");
    os.print(x2);
    os.print(" ");
    os.print(y2);
    os.print(" ");
    os.print(dim);
    os.println(" arcto");
    os.println("4 {pop} repeat");    
  }


  /** 
   * Draws an outlined rounded corner rectangle using the current color.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @param arcWidth the diameter of the arc
   * @param arcHeight the radius of the arc
   * @see #fillRoundRect
   */
  public void drawRoundRect(int x, int y, int width, int height, int arcWidth, int arcHeight) 
  {
    if (DEBUG) diagnostic("drawRoundRect(" + x + ", " +y + ", " +width + ", " +height + ", " + arcWidth + ", " + arcHeight +")");
    doRoundRect(x,y,width,height,arcWidth,arcHeight, false);
  }

  /** 
   * Draws a rounded rectangle filled in with the current color.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @param arcWidth the diameter of the arc
   * @param arcHeight the radius of the arc
   * @see #drawRoundRect
   */
  public void fillRoundRect(int x, int y, int width, int height, int arcWidth, int arcHeight) 
  {
    if (DEBUG) diagnostic("fillRoundRect(" + x + ", " +y + ", " +width + ", " +height + ", " +arcWidth + ", " +arcHeight+")");
    doRoundRect(x,y,width,height,arcWidth,arcHeight, true);
  }

  /**
   * Draws a highlighted 3-D rectangle.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @param raised a boolean that states whether the rectangle is raised or not
   */
  public void draw3DRect(int x, int y, int width, int height, boolean raised) 
  {
    if (DEBUG) diagnostic("draw3DRect(" + x + ", " +y + ", " +width + ", " +height + ", " +raised+")");
    Color c = getColor();
    Color brighter = c.brighter();
    Color darker = c.darker();

    setColor(raised ? brighter : darker);
    drawLine(x, y, x, y + height);
    drawLine(x + 1, y, x + width - 1, y);
    setColor(raised ? darker : brighter);
    drawLine(x + 1, y + height, x + width, y + height);
    drawLine(x + width, y, x + width, y + height);
    setColor(c);
  }    

  /**
   * Paints a highlighted 3-D rectangle using the current color.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @param raised a boolean that states whether the rectangle is raised or not
   */
  public void fill3DRect(int x, int y, int width, int height, boolean raised) 
  {
    if (DEBUG) diagnostic("fill3DRect(" + x + ", " +y + ", " +width + ", " +height + ", " +raised+")");
    Color c = getColor();
    Color brighter = c.brighter();
    Color darker = c.darker();

    if (!raised) 
      {
        setColor(darker);
      }
    fillRect(x+1, y+1, width-2, height-2);
    setColor(raised ? brighter : darker);
    drawLine(x, y, x, y + height - 1);
    drawLine(x + 1, y, x + width - 2, y);
    setColor(raised ? darker : brighter);
    drawLine(x + 1, y + height - 1, x + width - 1, y + height - 1);
    drawLine(x + width - 1, y, x + width - 1, y + height - 1);
    setColor(c);
  }    

  /** 
   * Draws an oval inside the specified rectangle using the current color.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @see #fillOval
   */
  public void drawOval(int x, int y, int width, int height) 
  {
    if (DEBUG) diagnostic("drawOval(" + x + ", " +y + ", " +width + ", " +height +")");
    doArc(x,y,width,height,0,360,false);
  }

  /** 
   * Fills an oval inside the specified rectangle using the current color.
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @see #drawOval
   */
  public void fillOval(int x, int y, int width, int height) 
  {
    if (DEBUG) diagnostic("fillOval(" + x + ", " +y + ", " +width + ", " +height+")");

    if(width == height) doArc(x,y,width);
    else                doArc(x,y,width,height,0,360,true);
    
  }

 
  protected void doArc(int x, int y, int width, int height,
                     int startAngle, int arcAngle, boolean fill)
  {
    if (DEBUG) diagnostic("doArc(" + x + ", " +y + ", " +width + ", " +height +
                          startAngle + ", " +arcAngle + ", " +fill+")");
    y = transformY(y);
    gsave();

    // cx,cy is the center of the arc
    int cx = x + width/2;
    int cy = y - height/2;

    // translate the page to be centered there
    translate(cx, cy);
    
    // scale the coordinate system - this is the only way to directly draw
    // an eliptical arc in postscript. Calculate the scale:
    
    float yscale = (float) height/(float)width;
    scale(1.0, yscale);
    
    if (fill) 
      moveto(0, 0);

    // now draw the arc.
    os.println("0 sl ");
    float endAngle = startAngle + arcAngle;
    os.print("0 0 ");
    os.print((float)width/2.0);
    os.print(" ");
    os.print(startAngle);
    os.print(" ");
    os.print(endAngle);
    os.println(" arc");

    if (fill)       
//      os.println("closepath");
      os.print(" C ");
    
    stroke(fill);
//    os.println("1 setlinewidth");
    os.print("1 sl ");
    
    // undo all the scaling!
    grestore();

  }
 

  protected void doArc(int x, int y, int width)
  {
    y = transformY(y);
    gsave();

    // cx,cy is the center of the arc
    int cx = x + width/2;
    int cy = y - width/2;

    // translate the page to be centered there
    translate(cx, cy);
    
    // scale the coordinate system - this is the only way to directly draw
    // an eliptical arc in postscript. Calculate the scale:
    
//    moveto(0, 0);

    // now draw the arc.
//    os.println("0 sl ");
//    os.print("0 0 ");
    os.print(" A ");
    os.print((float)width/2.0);
    os.println(" AA ");
//    os.print(" ");
//    os.print(0);
//    os.print(" ");
//    os.print(360);
//    os.println(" arc");

//    os.print(" C ");
    
//    stroke(true);

//    os.print("1 sl ");
    
    // undo all the scaling!
//    grestore();

  }
  
  
  /**
   * Draws an arc bounded by the specified rectangle from startAngle to
   * endAngle. 0 degrees is at the 3-o'clock position.Positive arc
   * angles indicate counter-clockwise rotations, negative arc angles are
   * drawn clockwise. 
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @param startAngle the beginning angle
   * @param arcAngle the angle of the arc (relative to startAngle).
   * @see #fillArc
   */
  public void drawArc(int x, int y, int width, int height,
                      int startAngle, int arcAngle) 
  {
    if (DEBUG) diagnostic("drawArc(" + x + ", " +y + ", " +width + ", " +height + ", " +startAngle + ", " +arcAngle+")");
    doArc(x,y,width,height,startAngle,arcAngle,false);
  }

  /** 
   * Fills an arc using the current color. This generates a pie shape.
   *
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the arc
   * @param height the height of the arc
   * @param startAngle the beginning angle
   * @param arcAngle the angle of the arc (relative to startAngle).
   * @see #drawArc
   */
  public void fillArc(int x, int y, int width, int height, int startAngle, int arcAngle) 
  {
    if (DEBUG) diagnostic("fillArc(" + x + ", " +y + ", " +width + ", " +height + ", " +startAngle + ", " +arcAngle+")");
    doArc(x,y,width,height,startAngle,arcAngle,true);
  }


  protected void doPoly(int xPoints[], int yPoints[], int nPoints, boolean fill, boolean close) 
  {
    if (DEBUG) diagnostic("doPoly(" + xPoints.length + ", " +yPoints.length + ", " +nPoints + ", " +fill + ", " +close+")");

    if (nPoints < 2)
      return;

    int newYPoints[] = new int[nPoints];
    int i;

    for (i=0; i< nPoints; i++)
      newYPoints[i] = transformY(yPoints[i]);

    moveto(xPoints[0], newYPoints[0]);

    for (i=0; i<nPoints; i++) 
      lineto(xPoints[i], newYPoints[i]);
    
    stroke(fill);
  }


  /** 
   * Draws a polygon defined by an array of x points and y points.
   * @param xPoints an array of x points
   * @param yPoints an array of y points
   * @param nPoints the total number of points
   * @see #fillPolygon
   */
  public void drawPolyline(int xPoints[], int yPoints[], int nPoints) 
  {
    if (DEBUG) diagnostic("drawPolyline(" + xPoints.length + ", " +yPoints.length + ", " +nPoints+")");
    doPoly(xPoints, yPoints, nPoints, false, false);
  }

  
  /** 
   * Draws a polygon defined by an array of x points and y points.
   * @param xPoints an array of x points
   * @param yPoints an array of y points
   * @param nPoints the total number of points
   * @see #fillPolygon
   */
  public void drawPolygon(int xPoints[], int yPoints[], int nPoints) 
  {
    if (DEBUG) diagnostic("drawPolygon(" + xPoints.length + ", " +yPoints.length + ", " +nPoints+")");
    doPoly(xPoints, yPoints, nPoints, false, true);
  }

  /** 
   * Draws a polygon defined by the specified point.
   * @param p the specified polygon
   * @see #fillPolygon
   */
  public void drawPolygon(Polygon p) 
  {
    if (DEBUG) diagnostic("drawPolygon(" + p + ")");
    doPoly(p.xpoints, p.ypoints, p.npoints, false, true);
  }
  
  /** 
   * Fills a polygon with the current color.
   * @param xPoints an array of x points
   * @param yPoints an array of y points
   * @param nPoints the total number of points
   * @see #drawPolygon
   */
  public void fillPolygon(int xPoints[], int yPoints[], int nPoints) 
  {
    if (DEBUG) diagnostic("fillPolygon(" + xPoints.length + ", " +yPoints.length + ", " +nPoints+")");
    doPoly(xPoints, yPoints, nPoints, true, true);
  }

  /** 
   * Fills the specified polygon with the current color.
   * @param p the polygon
   * @see #drawPolygon
   */
  public void fillPolygon(Polygon p) 
  {
    if (DEBUG) diagnostic("fillPolygon(" + p + ")");
    doPoly(p.xpoints, p.ypoints, p.npoints, true, true);
  }

  /** 
   * Draws the specified String using the current font and color.
   * The x,y position is the starting point of the baseline of the String.
   * @param str the String to be drawn
   * @param x the x coordinate
   * @param y the y coordinate
   * @see #drawChars
   * @see #drawBytes
   */

  public void drawString(String str, int x, int y) 
  {
    if (DEBUG) diagnostic("drawString(" + str + ", " +x + ", " +y +")");
    y = transformY(y);
    moveto(x, y);
    os.print(" (");
    os.print(str);
    os.println(") center show stroke");
  }


  /** 
   * Draws the specified characters using the current font and color.
   * @param data the array of characters to be drawn
   * @param offset the start offset in the data
   * @param length the number of characters to be drawn
   * @param x the x coordinate
   * @param y the y coordinate
   * @see #drawString
   * @see #drawBytes
   */
  public void drawChars(char data[], int offset, int length, int x, int y) 
  {
    if (DEBUG) diagnostic("drawChars(" + data.length + ", " +offset + ", " +length + ", " +x + ", " +y +")");
    drawString(new String(data, offset, length), x, y);
  }

  /** 
   * Draws the specified bytes using the current font and color.
   * @param data the data to be drawn
   * @param offset the start offset in the data
   * @param length the number of bytes that are drawn
   * @param x the x coordinate
   * @param y the y coordinate
   * @see #drawString
   * @see #drawChars
   */
  public void drawBytes(byte data[], int offset, int length, int x, int y) 
  {
    if (DEBUG) diagnostic("drawBytes(" + data.length + ", " +offset + ", " +length + ", " +x + ", " +y +")");
    drawString(new String(data, offset, length), x, y);
  }


  protected boolean doImage(Image img, int x, int y, int width, int height,
                          int sx, int sy, int sw, int sh,
                          ImageObserver observer, Color bgcolor) 
  {
    if (DEBUG) diagnostic("doImage(" + img + ", " +x + ", " +y + ", " +width + ", " + height + ", " + sx + ", " + sy + ", " + sw + ", " + sh + ", " + observer + ", " +bgcolor +")");
    
    int imgWidth = img.getWidth(observer);
    int imgHeight = img.getHeight(observer);

    y = transformY(y);

    int[] pix = new int[ imgWidth * imgHeight ];
    PixelGrabber pg = new PixelGrabber(img, 0, 0, imgWidth, imgHeight, pix, 0, imgWidth);

    boolean result = false;
    try
      {      
        result = pg.grabPixels();
      }
    catch (InterruptedException ie)
      {
        // FALL THROUGH
      }
    finally
      {
        if (!result)
          {
            os.println("%warning: error on image grab");
            System.err.println("warning: error on image grab: " + pg.getStatus());
            return false;
          }
      }


    // compute image size. First of all, if width or height is 0, image is 1:1.
    if (height < 1 || width < 1) 
      {
        height = imgHeight;
        width = imgWidth;
      }       

    int iLower = (sy == 0) ? 0 : sx;
    int iUpper = (sh == 0) ? imgHeight : sy + sh;

    int jLower = (sx == 0) ? 0 : sx;
    int jUpper = (sw == 0) ? imgWidth : sx + sw;
    
    int numYPixels = iUpper - iLower;
    int numXPixels = jUpper - jLower;

    gsave();

    os.println("% build a temporary dictionary");
    os.println("20 dict begin");
    emitColorImageProlog(numXPixels);

    os.println("% lower left corner");
    translate(x, y);

    os.println("% size of image");
    scale(width, height);

    os.print(numXPixels);
    os.print(" ");
    os.print(numYPixels);
    os.println(" 8");

    os.print("[");
    os.print(numXPixels);
    os.print(" 0 0 -");
    os.print(numYPixels);
    os.print(" 0 ");
    os.print(0);
    os.println("]");

    os.println("{currentfile pix readhexstring pop}");
    os.println("false 3 colorimage");
    os.println("");


    int offset, sleepyet=0;;
    // array to hold a line of pixel data
    char[] sb = new char[charsPerRow + 1];

    int bg = (bgcolor == null) ? -1 : bgcolor.getRGB();

    
    for (int i=iLower; i<iUpper; i++) 
      {
        offset = 0;
        ++sleepyet;

        // real color image. We're deliberately duplicating code here
        // in the interest of speed - we don't want to check bgcolor
        // on every iteration.
        for (int j=jLower; j<jUpper; j++) 
          {
            int coord = i * imgWidth + j;
            int n = pix[coord];
                        
            int alpha = n & 0xFF000000;

            if (alpha == 0)
              n = bg;

            // put hex chars into string
            // flip red for blue, to make postscript happy.

            sb[offset++] = hd[(n & 0xF00000) >> 20];
            sb[offset++] = hd[(n & 0xF0000)  >> 16];
            sb[offset++] = hd[(n & 0xF000)   >> 12];
            sb[offset++] = hd[(n & 0xF00)    >>  8];
            sb[offset++] = hd[(n & 0xF0)     >>  4];
            sb[offset++] = hd[(n & 0xF)           ];

            if (offset >= charsPerRow) 
              {
                os.write(sb, 0, offset);
                os.println();
                if (sleepyet > 5) 
                  {
                    try 
                      {
                        // let the screen update occasionally!
                        Thread.sleep(5);
                      } 
                    catch (InterruptedException ex) 
                      {
                        // yeah, so?
                      }
                    sleepyet = 0;
                  }
                offset = 0;
              }
          }
        
        // print partial rows
        if (offset != 0) 
          {            
            os.write(sb, 0, offset);
            os.println();
          }
      }
    
    os.println();
    os.println("end");
    grestore();
    
    return true;
  }
  
  /** 
   * Draws the specified image at the specified coordinate (x, y). If the image is 
   * incomplete the image observer will be notified later.
   * @param img the specified image to be drawn
   * @param x the x coordinate
   * @param y the y coordinate
   * @param observer notifies if the image is complete or not
   * @see Image
   * @see ImageObserver
   */

  public boolean drawImage(Image img, int x, int y,
                           ImageObserver observer) 
  {
    if (DEBUG) diagnostic("drawImage(" + img + ", " +x + ", " +y + ", " +observer+")");

    return doImage(img, x, y, 0, 0, 0, 0, 0, 0, observer, null);

  }
  
  /**
   * Warning this is not yet supported
   */
  public boolean drawImage(Image img, int x1, int y1,
                           int x2, int y2, 
                           int x3, int y3, int x4, int y4, 
                           ImageObserver observer) 
  {
    if (DEBUG) diagnostic("drawImage(" + img + ", " +x1 + ", " +y1 + ", " +x2 + ", " +y2 + ", " +x3 + ", " +y3 + ", " +x4 + ", " +y4 + ", " +observer +")"); 
    return doImage(img, x1, y1, x2-x1, y2-y1, x3, y3, x4-x3, y4-y3, observer, null);
  }

  /**
   * Warning this is not yet supported
   */
  public boolean drawImage(Image img, int x1, int y1,
                           int x2, int y2, 
                           int x3, int y3, int x4, int y4, 
                           Color c, 
                           ImageObserver observer) 
  {
    if (DEBUG) diagnostic("drawImage(" + img + ", " +x1 + ", " +y1 + ", " +x2 + ", " +y2 + ", " +x3 + ", " +y3 + ", " +x4 + ", " +y4 + ", " +c + ", " +observer+")");
    return doImage(img, x1, y1, x2-x1, y2-y1, x3, y3, x4-x3, y4-y3, observer, c);
  }

  /**
   * Draws the specified image inside the specified rectangle. The image is
   * scaled if necessary. If the image is incomplete the image observer will be
   * notified later.
   * @param img the specified image to be drawn
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @param observer notifies if the image is complete or not
   * @see Image
   * @see ImageObserver
   */
  public boolean drawImage(Image img, int x, int y,
                           int width, int height, 
                           ImageObserver observer) 
  {
    if (DEBUG) diagnostic("drawImage(" + img +"," +x + ", " +y + ", " +width + ", " +height + ", " + observer +")");
    return doImage(img, x, y, width, height, 0, 0, 0, 0, observer, null);
  }

  /** 
   * Draws the specified image at the specified coordinate (x, y). If the image is 
   * incomplete the image observer will be notified later.
   * @param img the specified image to be drawn
   * @param x the x coordinate
   * @param y the y coordinate
   * @param bgcolor the background color
   * @param observer notifies if the image is complete or not
   * @see Image
   * @see ImageObserver
   */

  public boolean drawImage(Image img, int x, int y, Color bgcolor,
                           ImageObserver observer) 
  {
    if (DEBUG) diagnostic("drawImage("+ img + ", " +x + ", " +y + ", " +bgcolor + ", " + observer +")");
    return doImage(img, x, y, 0, 0, 0, 0, 0, 0, observer, bgcolor);
  }

  /**
   * Draws the specified image inside the specified rectangle. The image is
   * scaled if necessary. If the image is incomplete the image observer will be
   * notified later.
   * @param img the specified image to be drawn
   * @param x the x coordinate
   * @param y the y coordinate
   * @param width the width of the rectangle
   * @param height the height of the rectangle
   * @param bgcolor the background color
   * @param observer notifies if the image is complete or not
   * @see Image
   * @see ImageObserver
   * NOTE: PSGr ignores the background color.
   */
  public boolean drawImage(Image img, int x, int y,
                           int width, int height, Color bgcolor,
                           ImageObserver observer) 
  {
    if (DEBUG) diagnostic("drawImage(" + img + ", " +x + ", " +y + ", " +width + ", " +height + ", " +bgcolor + ", " +observer +")");
    return doImage(img, x, y, width, height, 0, 0, 0, 0, observer, bgcolor);
  }
  
  /**
   * Disposes of this graphics context.  The Graphics context cannot be used after 
   * being disposed of.
   * @see #finalize
   */
  public void dispose() 
  {
    if (DEBUG) diagnostic("dispose() ");
    os.flush();
  }

  /**
   * Disposes of this graphics context once it is no longer referenced.
   * @see #dispose
   */
  public void finalize() 
  {
    if (DEBUG) diagnostic("finalize() ");
  }

  /**
   * Returns a String object representing this Graphic's value.
   */
  public String toString() 
  {
    if (DEBUG) diagnostic("toString() ");
    return getClass().getName() + "[font=" + getFont() + ",color=" + getColor() + "]";
  }

  /**
     Flip Y coords so Postscript looks like Java
  */

  protected int transformY(int y) 
  {
    return PAGEHEIGHT - y;
  }

  protected double transformY(double y) 
  {
    return PAGEHEIGHT - y;
  }
 
  
  /**
     Top of every PS file
  */

  protected void emitProlog() 
  {
//    os.println("%!PS-Adobe-2.0 EPSF-1.2 Created by PSGr Java PostScript Context");
//    os.println("%%BoundingBox:0 0 595 840");
//    os.println("%!PS-Adobe-3.0");      
//    os.println("%%BoundingBox:30 30 565 810");    
    os.println("% PSGr is (C) 1996 Ernest Friedman-Hill and Sandia National Labs");
    os.println("% Right to unrestricted personal and commerical use is granted");
    os.println("% if this acknowledgement is given on product or packing materials");
    os.println("/bd {bind def} def");
    os.println("/M {moveto} bd");
    os.println("/L {lineto} bd");
    os.println("/S {stroke} bd");
    os.println("/C {closepath} bd");
    os.println("/T {translate} bd");
    os.println("/sc {scale} bd");
    os.println("/gr {grestore} bd");
    os.println("/gs {gsave} bd");
    os.println("/sl {setlinewidth} bd");
    os.println("/E {eofill} bd");
    os.println("/A {0 0 moveto 0 sl 0 0} bd");
    os.println("/AA {0 360 arc closepath gsave eofill grestore stroke 1 setlinewidth grestore} bd");
    os.println("/center {dup stringwidth pop 2 div neg 0 rmoveto} bd");
    
    translate(XOFFSET, -YOFFSET);
    setFont(font);
  }


  protected void emitColorImageProlog(int xdim) 
  {
    os.println("% Color picture stuff, lifted from XV's PS files");

    os.println("% define string to hold a scanline's worth of data");
    os.print("/pix ");
    os.print(xdim*3);
    os.println(" string def");

    os.println("% define space for color conversions");
    os.print("/grays ");
    os.print(xdim);
    os.println(" string def  % space for gray scale line");
    os.println("/npixls 0 def");
    os.println("/rgbindx 0 def");

    os.println("% define 'colorimage' if it isn't defined");
    os.println("%   ('colortogray' and 'mergeprocs' come from xwd2ps");
    os.println("%     via xgrab)");
    os.println("/colorimage where   % do we know about 'colorimage'?");
    os.println("{ pop }           % yes: pop off the 'dict' returned");
    os.println("{                 % no:  define one");
    os.println("/colortogray {  % define an RGB->I function");
    os.println("/rgbdata exch store    % call input 'rgbdata'");
    os.println("rgbdata length 3 idiv");
    os.println("/npixls exch store");
    os.println("/rgbindx 0 store");
    os.println("0 1 npixls 1 sub {");
    os.println("grays exch");
    os.println("rgbdata rgbindx       get 20 mul    % Red");
    os.println("rgbdata rgbindx 1 add get 32 mul    % Green");
    os.println("rgbdata rgbindx 2 add get 12 mul    % Blue");
    os.println("add add 64 idiv      % I = .5G + .31R + .18B");
    os.println("put");
    os.println("/rgbindx rgbindx 3 add store");
    os.println("} for");
    os.println("grays 0 npixls getinterval");
    os.println("} bind def");
    os.println("");
    os.println("% Utility procedure for colorimage operator.");
    os.println("% This procedure takes two procedures off the");
    os.println("% stack and merges them into a single procedure.");
    os.println("");
    os.println("/mergeprocs { % def");
    os.println("dup length");
    os.println("3 -1 roll");
    os.println("dup");
    os.println("length");
    os.println("dup");
    os.println("5 1 roll");
    os.println("3 -1 roll");
    os.println("add");
    os.println("array cvx");
    os.println("dup");
    os.println("3 -1 roll");
    os.println("0 exch");
    os.println("putinterval");
    os.println("dup");
    os.println("4 2 roll");
    os.println("putinterval");
    os.println("} bind def");
    os.println("");
    os.println("/colorimage { % def");
    os.println("pop pop     % remove 'false 3' operands");
    os.println("{colortogray} mergeprocs");
    os.println("image");
    os.println("} bind def");
    os.println("} ifelse          % end of 'false' case");

  }

  public void gsave() 
  {
//    os.println("gsave");
      os.print(" gs ");
  }

  public void grestore() 
  {
//    os.println("grestore");
    os.print(" gr ");      
  }

  public void emitThis(String s) 
  {
    os.println(s);
  }

  protected void diagnostic(String s)
  {
    os.print("% PSGR-");
    os.print(hashCode());
    os.print(": ");
    os.println(s);
  }
  
}


  



