import java.io.File;
import javax.swing.*;
import javax.swing.filechooser.*;

/**
 *OutputFilter2.java
 *
 *Created on 06.02.2002
 *
 *Last modified on 06.02.2002
 */

public class OutputFilter2 extends FileFilter {
    
public boolean accept(File f)
	{
	if (f.isDirectory())
		{
		return true;
		}

	String extension = getExtension(f);
	if (extension != null)
		{
		if (extension.equals("txt"))
			{
  		                return true;
     			}
		else
			{
		                return false;
			}
		}
	return false;
    	}

public static String getExtension(File f)
	{
	String ext = null;
	String s = f.getName();
	int i = s.lastIndexOf('.');
	if (i > 0 &&  i < s.length() - 1) 
		{
		ext = s.substring(i+1).toLowerCase();
		}
	return ext;
	}

public String getDescription()
	{
	return "Text File (*.txt)";
	}

}
