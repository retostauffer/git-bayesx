//package gov.sandia.postscript;

import java.awt.Graphics;
import java.io.Writer;
import java.text.AttributedCharacterIterator;

/**
 * PSGr2 is a Graphics subclass for Java 2 that images to PostScript.
 * (C) 1996 E.J. Friedman-Hill and Sandia National Labs
 * @version 	2.1
 * @author 	Ernest Friedman-Hill
 * @author      ejfried@ca.sandia.gov
 * @author      http://herzberg.ca.sandia.gov
 */

public class PSGr2 extends PSGrBase
{
  /**
   * Constructs a new PSGr2 Object. Unlike regular Graphics objects,
   * PSGr contexts can be created directly.
   * @param o Output stream for PostScript output
   * @see #create
   */

  public PSGr2()
  {
    super();
  }

  /**
   * Constructs a new PSGr2 Object. Unlike regular Graphics objects,
   * PSGr contexts can be created directly.
   * @param o Output stream for PostScript output
   * @see #create
   */

  public PSGr2(Writer o)
  {
    super(o, true);
  }

  /**
   * Constructs a new PSGr2 Object. Unlike regular Graphics objects,
   * PSGr contexts can be created directly.
   * @param o Output stream for PostScript output
   * @see #create
   */

  public PSGr2(Writer o, Graphics g, boolean emitProlog)
  {
    super(o, emitProlog);
  }

  /**
   * So far Unimplemented Java 2 addition.
   */
  public void drawString(AttributedCharacterIterator i, int x, int y)
  {
    throw new RuntimeException("drawString(AttributedCharacterIterator, int, int) not implemented");
  }

}


  
