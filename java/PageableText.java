import java.awt.*;
import java.awt.print.*;
import java.io.*;
import java.util.Vector;

public class PageableText implements Pageable, Printable
{

public static String FONTFAMILY = "Monospaced";
public static int FONTSIZE = 10;
public static int FONTSTYLE = Font.PLAIN;
public static float LINESPACEFACTOR = 1.1f;

PageFormat format;
Vector lines;
Font font;
int linespacing;
int linesPerPage;
int numPages;
int baseline = -1;

public PageableText(File file, PageFormat format) throws IOException
	{
	this(new FileReader(file), format);
	}

public PageableText(Reader stream, PageFormat format) throws IOException
	{
	this.format = format;
	BufferedReader in = new BufferedReader(stream);
	lines = new Vector();
	String line;
	while((line=in.readLine()) != null)
		{
		lines.addElement(line);
		}
	font = new Font(FONTFAMILY, FONTSTYLE, FONTSIZE);
	linespacing = (int)(FONTSIZE*LINESPACEFACTOR);
	linesPerPage = (int)Math.floor(format.getImageableHeight()/linespacing);
	numPages = (lines.size()-1)/linesPerPage + 1;
	}
	

public int getNumberOfPages()
	{
	return numPages;
	}

public PageFormat getPageFormat(int pagenum)
	{
	return format;
	}

public Printable getPrintable(int pagenum)
	{
	return this;
	}

public int print(Graphics g, PageFormat format, int pagenum)
	{
	if((pagenum<0) | (pagenum>=numPages))
		{
		return NO_SUCH_PAGE;
		}
	if(baseline==-1)
		{
		FontMetrics fm = g.getFontMetrics(font);
		baseline = fm.getAscent();
		}
	g.setColor(Color.white);
	g.fillRect((int)format.getImageableX(), (int)format.getImageableY(), (int)format.getImageableWidth(), (int)format.getImageableHeight());
	g.setFont(font);
	g.setColor(Color.black);
	int startLine = pagenum * linesPerPage;
	int endLine = startLine + linesPerPage-1;
	if(endLine >= lines.size())
		{
		endLine = lines.size()-1;
		}
	int x0 = (int)format.getImageableX();
	int y0 = (int)format.getImageableY()+baseline;
	for(int i=startLine; i<=endLine; i++)
		{
		String line = (String)lines.elementAt(i);
		if(line.length()>0)
			{
			g.drawString(line, x0, y0);
			}
		y0 += linespacing;
		}
	return PAGE_EXISTS;
	}
}
