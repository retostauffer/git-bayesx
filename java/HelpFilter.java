import java.io.File;
import javax.swing.*;
import javax.swing.filechooser.*;

/**
 *OutputFilter.java
 *
 *Created on 26.02.2001
 *
 *Last modified on 26.02.2001
 */

public class HelpFilter extends FileFilter {
    
public boolean accept(File f)
	{
	if (f.isDirectory())
		{
		return true;
		}

	String extension = getExtension(f);
	if (extension != null)
		{
		if (extension.equals("exe"))
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
	return "Executable (*.exe)";
	}

}
