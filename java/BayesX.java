import javax.swing.*;
import javax.swing.text.*;
import javax.swing.text.rtf.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.print.*;
import java.awt.datatransfer.*;
import java.util.*;
import java.io.*;
import java.beans.*;

import umontreal.iro.lecuyer.rng.RandMrg;
import umontreal.iro.lecuyer.randvar.UniformGen;
import umontreal.iro.lecuyer.randvar.NormalGen;
import umontreal.iro.lecuyer.randvar.ExponentialGen;
import umontreal.iro.lecuyer.randvar.BinomialGen;
import umontreal.iro.lecuyer.randvar.GammaGen;
import umontreal.iro.lecuyer.randvar.PoissonGen;
import umontreal.iro.lecuyer.randvar.WeibullGen;
//import umontreal.iro.lecuyer.randvar.*;
//import umontreal.iro.lecuyer.probdist.*;

//import gov.sandia.postscript.PSGrBase;
//import gov.sandia.postscript.PSGr1;
//import gov.sandia.postscript.PSGr2;

//import javax.help.*;
import java.net.*;

/**
 *BayesX.java
 *
 *Created on 14.02.2001
 *
 *Last modified on 17.08.2004
 *
 *@author Thomas Kneib
 */

public class BayesX extends JFrame
	implements ActionListener, ItemListener, DocumentListener, ListSelectionListener,
		WindowListener, AdjustmentListener//, CaretListener
{

//Funktionen, die die Kommunikation mit C++ steuern

static
	{
	System.loadLibrary("BayesXdll");
	}


private String defaultDirectory = (new File(System.getProperty("user.dir"),"output")).toString();

//Variablen der Klasse BayesX
private JMenuBar jMenuBar;
	private JMenu file;
		private final JFileChooser fileChooser;
		private JMenuItem clear;
		private JMenuItem open;
		private JMenuItem save;
		private JMenuItem saveas;
		private JMenuItem print;
		private JMenuItem exit;
	private JMenu preferences;
		private ButtonGroup group2;
			private JCheckBoxMenuItem small;
			private JCheckBoxMenuItem normal;
			private JCheckBoxMenuItem large;
				short scriptsize;
				double fontFactor;
			private JMenuItem defwin;
/*		private JCheckBoxMenuItem bold;
		private JCheckBoxMenuItem italic;
		private JCheckBoxMenuItem underline;*/
	private JMenu window;
		private ButtonGroup group;
			private JCheckBoxMenuItem comm;
			private JCheckBoxMenuItem out;
			private JCheckBoxMenuItem rev;
			private JCheckBoxMenuItem obj;
	private JMenu help;
	private JMenuItem help1;
	private JMenuItem help2;
	private JMenuItem help3;
	private JMenuItem info;
/*	private JMenuItem help2;
	private File helpReader;
	private String helpString;
	private JFileChooser fileChooser1;*/

private JDesktopPane jDesktopPanel;
	private JInternalFrame output;
		private boolean hasBeenSaved;
		private boolean isSaved;
		private JEditorPane outputPane;
		private RTFEditorKit outputEditor;
		private File outputFile;
		private DefaultStyledDocument outputDocument;
		private SimpleAttributeSet outputAttributes;
		private JScrollPane outputScrollPane;
		private File printFile;
		private DefaultEditorKit printEditor;
	private JInternalFrame command;
		private char delimiter;
		private JTextArea commandArea;
		private JScrollPane commandScrollPane;
		private ParseThread parseThread = new ParseThread(this);
		private int threadPriority = Thread.NORM_PRIORITY;
	private JInternalFrame review;
		private JList reviewList;
			private Vector reviewVector;
		private JScrollPane reviewScrollPane;
	private JInternalFrame objects;
		private JSplitPane objectsSplitPane;
		private JList objectsList1;
			private Vector objectsVector;
		private JList objectsList2;
		private JScrollPane objectsScrollPane;
			private Vector objectVector;
		private JFrame objectFrame;
			private Vector varnames;
			private int rows;
			private int cols;
			private Object[][] data;
			private Object[] vars;
			private Object[][] rowHeaderData1;
			private Object[] rowHeaderData2;
			private DefaultTableModel tableModel1;
			private DefaultTableModel tableModel2;
			private JTable rowHeaderTable;
			private JTable objectTable;
			private JScrollBar verticalScrollBar;
			private JScrollPane horizontalScrollPane;
		private JFrame mapFrame;
			protected MapPanel mapPanel;
			private JFileChooser fileChooser2;
			private JScrollPane mapFrameScrollPane;

private JPanel buttonPanel;
	private JButton breakButton;
	private JButton pauseButton;
	protected boolean pause;
	private boolean processRunning;
	private JButton outputButton;
	private JLabel priorityLabel;
	private JComboBox priorityBox;
private File registryFile;
	private int[] registryArray;
	private String registryString;
private boolean consoleInput;

// Plot-Optionen
protected short function;
protected short plotsperpage;

// Optionen fuer Javadrawmap
protected boolean color;
protected boolean legend;
protected boolean swap;
protected boolean drawnames;
protected boolean hcl;
protected String title;
protected String outfile;
protected double upperlimit;
protected double lowerlimit;
protected short shades;

protected int nrNA;

// Optionen fuer Javaplotnonp
protected String xlab;
protected String ylab;
protected String connect;
protected String linecolor;
protected int width;
protected int height;
protected double xmax;
protected double xmin;
protected double ymax;
protected double ymin;
protected double xstep;
protected double ystep;
protected double xstart;
protected double ystart;
protected int year;
protected int month;
protected int linewidth;
protected int fontsize;
protected int pointsize;
protected double titlescale;

// Optionen fuer Javaplotautocor
protected boolean meanautocor;

// Optionen fuer Javaplotsurf
protected String zlab;
protected double xrot;
protected double yrot;
protected double zrot;
protected double zmin;
protected double zmax;
protected double zstart;
protected double zstep;
protected int gridsize;

// Random Number Generation

RandMrg rStream;

//Konstruktor zur Erzeugung eines BayesX-Fensters
public BayesX()
	{

//Setze consoleInput per default auf false => Version mit Fenstern in JavaOutput

	consoleInput = false;

//Lese bzw. erzeuge die Datei, die die Groessen der Fenster enthaelt

	registryFile = new File(new File(System.getProperty("user.dir")),"registry.bayesx");
	try
		{
		if(!registryFile.exists() )
			{
			PrintWriter out = new PrintWriter(new FileWriter(registryFile));
			if(out!=null)
				{
				out.println("750 500 5 115 450 285 5 5 450 105 460 5 275 200 460 210 275 190 12\nBayesX System-File. DO NOT EDIT!");
				out.close();
				}
			}
		else
			{
//			System.out.println("write check failed 1");
			}
		BufferedReader in = new BufferedReader(new FileReader(registryFile.getName()));
		if(registryFile.canWrite() && in!=null)
			{
			registryString = in.readLine();
			in.close();
			registryArray = new int[19];
			StringTokenizer st = new StringTokenizer(registryString);
			for(int i=0; i<19; i++)
				{
				Integer k = new Integer(st.nextToken());
				registryArray[i] = k.intValue();
				}
			}
		else
			{
//			System.out.println("write check failed 2");
			registryArray = new int[19];
			int[] regArray = {750,500,5,115,450,285,5,5,450,105,460,5,275,200,460,210,275,190,12};
			for(int i=0;i<regArray.length; i++)
				{
//				System.out.println(i);
				registryArray[i] = regArray[i];
				}
//			System.out.println("write check failed 3");
			}
		}
	catch(IOException ioe)
		{
		System.err.println(ioe.getMessage());
		}

	this.getContentPane().setLayout(new BorderLayout(5,5));

//Erzeugen der Menueleiste

	jMenuBar = new JMenuBar();
	jMenuBar.setRequestFocusEnabled(false);
	fileChooser = new JFileChooser(defaultDirectory);
	fileChooser.addChoosableFileFilter(new OutputFilter2());
	fileChooser.addChoosableFileFilter(new OutputFilter());
	file = new JMenu("File");
	file.setRequestFocusEnabled(false);
	file.setMnemonic(KeyEvent.VK_F);
	clear = new JMenuItem("Clear Output Window");
	clear.setRequestFocusEnabled(false);
	clear.setMnemonic(KeyEvent.VK_C);
	clear.addActionListener(this);
	open = new JMenuItem("Open Output File");
	open.setRequestFocusEnabled(false);
	open.setMnemonic(KeyEvent.VK_O);
	open.addActionListener(this);
	save = new JMenuItem("Save Output");
	save.setRequestFocusEnabled(false);
	save.setMnemonic(KeyEvent.VK_S);
	save.addActionListener(this);
	save.setEnabled(false);
	saveas = new JMenuItem("Save Output as...");
	saveas.setRequestFocusEnabled(false);
	saveas.setMnemonic(KeyEvent.VK_A);
	saveas.addActionListener(this);
	print = new JMenuItem("Print Output");
	print.setRequestFocusEnabled(false);
	print.setMnemonic(KeyEvent.VK_P);
	print.addActionListener(this);
	exit = new JMenuItem("Exit");
	exit.setRequestFocusEnabled(false);
	exit.setMnemonic(KeyEvent.VK_E);
	exit.addActionListener(this);
	file.add(clear);
	file.add(open);
	file.add(save);
	file.add(saveas);
	file.add(print);
	file.add(exit);

	preferences = new JMenu("Preferences");
	preferences.setRequestFocusEnabled(false);
	preferences.setMnemonic(KeyEvent.VK_P);
	group2 = new ButtonGroup();
	small = new JCheckBoxMenuItem("Fontsize Small");
	small.addActionListener(this);
	small.setMnemonic(KeyEvent.VK_S);
	normal = new JCheckBoxMenuItem("Fontsize Normal");
	normal.addActionListener(this);
	normal.setMnemonic(KeyEvent.VK_N);
	large = new JCheckBoxMenuItem("Fontsize Large");
	large.addActionListener(this);
	large.setMnemonic(KeyEvent.VK_L);
	group2.add(small);
	group2.add(normal);
	group2.add(large);
	scriptsize=(short)registryArray[18];
	if(scriptsize==10)
		{
		fontFactor=0.84;
		}
	else if(scriptsize==12)
		{
		fontFactor=1.0;
		}
	else
		{
		fontFactor=1.4;
		}

	defwin = new JMenuItem("Default Windowing");
	defwin.addActionListener(this);
	defwin.setMnemonic(KeyEvent.VK_D);

	preferences.add(small);
	preferences.add(normal);
	preferences.add(large);
	preferences.add(defwin);


/*	bold = new JCheckBoxMenuItem("Output Bold");
	bold.setRequestFocusEnabled(false);
	bold.setMnemonic(KeyEvent.VK_B);
	bold.addActionListener(new StyledEditorKit.BoldAction());
	italic = new JCheckBoxMenuItem("Output Italic");
	italic.setRequestFocusEnabled(false);
	italic.setMnemonic(KeyEvent.VK_I);
	italic.addActionListener(new StyledEditorKit.ItalicAction());
	underline = new JCheckBoxMenuItem("Output Underline");
	underline.setRequestFocusEnabled(false);
	underline.setMnemonic(KeyEvent.VK_U);
	underline.addActionListener(new StyledEditorKit.UnderlineAction());
	edit.add(bold);
	edit.add(italic);
	edit.add(underline);*/

	window = new JMenu("Window");
	window.setRequestFocusEnabled(false);
	window.setMnemonic(KeyEvent.VK_W);
	group = new ButtonGroup();
	comm = new JCheckBoxMenuItem("Command");
	comm.setRequestFocusEnabled(false);
	comm.setMnemonic(KeyEvent.VK_C);
	comm.setSelected(true);
	comm.addItemListener(this);
	out = new JCheckBoxMenuItem("Output");
	out.setRequestFocusEnabled(false);
	out.setMnemonic(KeyEvent.VK_O);
	out.addItemListener(this);
	obj = new JCheckBoxMenuItem("ObjectBrowser");
	obj.setRequestFocusEnabled(false);
	obj.setMnemonic(KeyEvent.VK_B);
	obj.addItemListener(this);
	rev = new JCheckBoxMenuItem("Review");
	rev.setRequestFocusEnabled(false);
	rev.setMnemonic(KeyEvent.VK_R);
	rev.addItemListener(this);
	group.add(comm);
	group.add(obj);
	group.add(out);
	group.add(rev);
	window.add(comm);
	window.add(obj);
	window.add(out);
	window.add(rev);

	help = new JMenu("Help");
/*	fileChooser1 = new JFileChooser(defaultDirectory);
	fileChooser1.addChoosableFileFilter(new HelpFilter());*/
        help1 = new JMenuItem("Reference Manual");
        help1.setMnemonic(KeyEvent.VK_R);
        help.add(help1);
        help1.addActionListener(this);
        help2 = new JMenuItem("Tutorials Manual");
        help2.setMnemonic(KeyEvent.VK_T);
        help.add(help2);
        help2.addActionListener(this);
        help3 = new JMenuItem("Methodology Manual");
        help3.setMnemonic(KeyEvent.VK_M);
        help.add(help3);
        help3.addActionListener(this);
	info = new JMenuItem("About BayesX");
	info.setMnemonic(KeyEvent.VK_A);
	help.add(info);
	info.addActionListener(this);
/*	help2 = new JMenuItem("Change Help-Path");
	help2.setMnemonic(KeyEvent.VK_C);
	help2.addActionListener(this);
	help.add(help2);*/
	help.setMnemonic(KeyEvent.VK_H);

	jMenuBar.add(file);
	jMenuBar.add(preferences);
	jMenuBar.add(window);
	jMenuBar.add(help);

	setJMenuBar(jMenuBar);


//Und jetzt das Panel fuer die Buttons sowie die Buttons selbst

	buttonPanel = new JPanel(new FlowLayout(FlowLayout.LEFT,10,5));

	breakButton = new JButton("BREAK");
	breakButton.setToolTipText("Stops the current process");
	breakButton.addActionListener(this);
	breakButton.setRequestFocusEnabled(false);
	breakButton.setEnabled(false);

	pauseButton = new JButton("PAUSE");
	pause = false;
	processRunning = false;
	pauseButton.setToolTipText("Interrupts the current process");
	pauseButton.addActionListener(this);
	pauseButton.setRequestFocusEnabled(false);

	outputButton = new JButton("SUPPRESS OUTPUT");
	outputButton.setToolTipText("Suppresses the documentation of output");
	outputButton.addActionListener(this);
	outputButton.setRequestFocusEnabled(false);

	breakButton.setPreferredSize(new Dimension(100, 30));
	pauseButton.setPreferredSize(new Dimension(100, 30));
	outputButton.setPreferredSize(new Dimension(160, 30));

	priorityLabel = new JLabel("PRIORITY:");
	priorityLabel.setPreferredSize(new Dimension(58, 30));

	String[] priorityStrings = {" VERY LOW", " LOW", " NORMAL", " HIGH", " VERY HIGH"};
	priorityBox = new JComboBox(priorityStrings);
	priorityBox.setSelectedIndex(2);
	priorityBox.setPreferredSize(new Dimension(100, 30));
	priorityBox.setToolTipText("Change the priority of BayesX");
	priorityBox.setRequestFocusEnabled(false);
	priorityBox.addActionListener(this);

	buttonPanel.add(breakButton);
	buttonPanel.add(pauseButton);
	buttonPanel.add(outputButton);
	buttonPanel.add(priorityLabel);
	buttonPanel.add(priorityBox);
	this.getContentPane().add(buttonPanel,BorderLayout.NORTH);
	buttonPanel.setBorder(BorderFactory.createLineBorder(Color.black));

//Erzeugen des Desktops fuer die einzelnen Fenster

	jDesktopPanel = new JDesktopPane();
	getContentPane().add(jDesktopPanel,BorderLayout.CENTER);
	jDesktopPanel.setDragMode(JDesktopPane.OUTLINE_DRAG_MODE);

//Die einzelnen Fenster erzeugen und hinzufuegen

//Zuerst das Output-Fenster

	output = new JInternalFrame("Output");
	output.setVisible(true);
	output.setMaximizable(true);
	output.setIconifiable(true);
	output.setResizable(true);
	outputPane = new JEditorPane();
	outputPane.addMouseListener(new MouseAdapter()
		{
		public void mouseClicked(MouseEvent me)
			{
			try
				{
				out.setSelected(true);
				output.setSelected(true);
				}
			catch(PropertyVetoException pve)
				{
				System.err.println(pve.getMessage());
				}
			}
		});
//	outputPane.addCaretListener(this);
	printFile = new File(new File(System.getProperty("user.dir")),"print.bayesx");
	printEditor = new DefaultEditorKit();
	outputEditor = new RTFEditorKit();
	outputPane.setEditorKit(outputEditor);
	outputDocument = new DefaultStyledDocument();
	outputPane.setDocument(outputDocument);
/*	final Clipboard cb = Toolkit.getDefaultToolkit().getSystemClipboard();

	outputPane.addKeyListener(new KeyAdapter()
		{
		public void keyReleased(KeyEvent e)
			{
			if(e.getModifiers()==InputEvent.CTRL_MASK && outputPane.getSelectedText()!=null)
				{
				if(e.getKeyCode()==KeyEvent.VK_X)
					{
					StringSelection sel = new StringSelection(outputPane.getSelectedText());
					cb.setContents(sel,sel);
					outputPane.replaceSelection("");
					}
				else if(e.getKeyCode()==KeyEvent.VK_C)
					{
					StringSelection sel = new StringSelection(outputPane.getSelectedText());
					cb.setContents(sel,sel);
					}
				}
			}
		});*/
	outputAttributes = new SimpleAttributeSet();
	StyleConstants.setBold(outputAttributes,false);
	StyleConstants.setItalic(outputAttributes,false);
	StyleConstants.setFontSize(outputAttributes,12);
	StyleConstants.setForeground(outputAttributes,Color.black);
	StyleConstants.setFontFamily(outputAttributes, "Courier New");
	outputDocument.addDocumentListener(this);
	outputScrollPane = new JScrollPane(outputPane);
	outputScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
	output.getContentPane().add(outputScrollPane,BorderLayout.CENTER);
	output.addInternalFrameListener(new InternalFrameAdapter()
		{
		public void internalFrameActivated(InternalFrameEvent ife)
			{
			outputPane.requestFocus();
			if(!processRunning)
				{
				outputPane.getCaret().setVisible(true);
				commandArea.getCaret().setVisible(false);
				}
			out.setSelected(true);
			}
		});
	jDesktopPanel.add(output);
	output.setBounds(registryArray[2],registryArray[3],registryArray[4],registryArray[5]);
	Out("BayesX - Software for Bayesian Inference in Structured Additive Regression Models\n\n",true,false,(short)11,0,0,0);
	Out("Version 2.1 (07.05.2012)\n\n");

//	Out("Note: When running time consuming computations it is useful to reduce the priority of BayesX in the Windows Task-Manager!\n\n");
	hasBeenSaved = false;
	isSaved = true;

//Testen, ob man Schreibrechte in den Defaultverzeichnissen hat

	if(!(new File(System.getProperty("user.dir"),"output")).canWrite() || !(new File(System.getProperty("user.dir"),"temp")).canWrite())
		{
		Out("WARNING: No permission to write to default directories.\n",true,true,(short)11,255,0,0);
		Out("         Specify a new default directory using the defaultpath command.\n"+
		    "         Type for example: defaultpath=c:\\temp",false,false,(short)11,0,0,0);
		}

//und dann das Command-Fenster

	delimiter = '\n';
	command = new JInternalFrame("Command");
	command.setVisible(true);
	command.setMaximizable(true);
	command.setIconifiable(true);
	command.setResizable(true);
	commandArea = new JTextArea();
	commandArea.addKeyListener(new KeyAdapter()
		{
		public void keyPressed(KeyEvent e)
			{
			if(e.getKeyChar()==delimiter)
				{
				doparse(commandArea.getText());
				}
			}
		});
	commandArea.addMouseListener(new MouseAdapter()
		{
		public void mouseClicked(MouseEvent me)
			{
			try
				{
				comm.setSelected(true);
				command.setSelected(true);
				}
			catch(PropertyVetoException pve)
				{
				System.err.println(pve.getMessage());
				}
			}
		});
	commandArea.setLineWrap(true);
	commandArea.setWrapStyleWord(true);
	commandArea.setFont(new Font("Monospaced",Font.PLAIN,12));
	commandScrollPane = new JScrollPane(commandArea);
	commandScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
	command.getContentPane().add(commandScrollPane,BorderLayout.CENTER);
	command.addInternalFrameListener(new InternalFrameAdapter()
		{
		public void internalFrameActivated(InternalFrameEvent ife)
			{
			commandArea.requestFocus();
			if(!processRunning)
				{
				commandArea.getCaret().setVisible(true);
				outputPane.getCaret().setVisible(false);
				}
			comm.setSelected(true);
			}
		});
	jDesktopPanel.add(command);
	command.setBounds(registryArray[6],registryArray[7],registryArray[8],registryArray[9]);

//und das Review-Fenster

	reviewVector = new Vector();
	review = new JInternalFrame("Review");
	review.setVisible(true);
	review.setMaximizable(true);
	review.setIconifiable(true);
	review.setResizable(true);
	reviewList = new JList(reviewVector);
	reviewList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
	reviewList.addMouseListener(new MouseAdapter()
		{
		public void mouseClicked(MouseEvent me)
			{
			try
				{
				rev.setSelected(true);
				review.setSelected(true);
				}
			catch(PropertyVetoException pve)
				{
				System.err.println(pve.getMessage());
				}
			if(me.getClickCount()>1)
				{
				int index = reviewList.getSelectedIndex();
				if(index>-1)
					{
					String selected = (String)reviewVector.elementAt(index);
					commandArea.setText(selected);
					try
						{
						command.setSelected(true);
						}
					catch(PropertyVetoException pve)
						{
						System.err.println(pve.getMessage());
						}
					}
				}
			}
		});
	reviewList.addKeyListener(new KeyAdapter()
		{
		public void keyPressed(KeyEvent e)
			{
			if(e.getKeyChar()=='\n')
				{
				int index = reviewList.getSelectedIndex();
				if(index>-1)
					{
					String selected = (String)reviewVector.elementAt(index);
					commandArea.setText(selected);
					try
						{
						command.setSelected(true);
						}
					catch(PropertyVetoException pve)
						{
						System.err.println(pve.getMessage());
						}
					}
				}
			}
		});
	review.addInternalFrameListener(new InternalFrameAdapter()
		{
		public void internalFrameActivated(InternalFrameEvent ife)
			{
			reviewList.requestFocus();
			rev.setSelected(true);
			}
		public void internalFrameDeactivated(InternalFrameEvent ife)
			{
			if(ife.getSource() == review)
				{
				reviewList.clearSelection();
				}
			}
		});
	reviewScrollPane = new JScrollPane(reviewList);
	review.getContentPane().add(reviewScrollPane,BorderLayout.CENTER);
	jDesktopPanel.add(review);
	review.setBounds(registryArray[10],registryArray[11],registryArray[12],registryArray[13]);

//und das Objects-Fenster

	objectsVector = new Vector();
	setObjectTypeList(objectsVector);
	objectVector = new Vector();
	objects = new JInternalFrame("ObjectBrowser");
	objects.setVisible(true);
	objects.setMaximizable(true);
	objects.setIconifiable(true);
	objects.setResizable(true);
	objectsList1 = new JList(objectsVector);
	objectsList1.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
	objectsList1.addMouseListener(new MouseAdapter()
		{
		public void mouseClicked(MouseEvent me)
			{
			try
				{
				obj.setSelected(true);
				objects.setSelected(true);
				}
			catch(PropertyVetoException pve)
				{
				System.err.println(pve.getMessage());
				}
			int index = objectsList1.getSelectedIndex();
			if(index>-1)
				{
				objectVector.clear();
				String selected = (String)objectsVector.elementAt(index);
				setObjectList(objectVector,selected);
				objectsList2.setListData(objectVector);
				}
			}
		});
	objectsList1.addListSelectionListener(this);
	objectsList2 = new JList();
	objectsList2.addMouseListener(new MouseAdapter()
		{
		public void mouseClicked(MouseEvent me)
			{
			try
				{
				obj.setSelected(true);
				objects.setSelected(true);
				}
			catch(PropertyVetoException pve)
				{
				System.err.println(pve.getMessage());
				}
			if(me.getClickCount()>=2)
				{
				int index = objectsList2.getSelectedIndex();
				if(index>-1)
					{
					String str = (String)objectVector.elementAt(index)+".describe";
					doparse(str);
					}
				}
			}
		});
	objectsList2.addKeyListener(new KeyAdapter()
		{
		public void keyReleased(KeyEvent ke)
			{
			if(ke.getKeyCode()==KeyEvent.VK_ENTER)
				{
				int index = objectsList2.getSelectedIndex();
				if(index >-1)
					{
					String str = (String)objectVector.elementAt(index)+".describe";
					doparse(str);
					}
				}
			}
		public void keyPressed(KeyEvent e)
			{
			if(e.getKeyCode()==KeyEvent.VK_DELETE)
				{
				int index = objectsList2.getSelectedIndex();
				if(index >-1)
					{
					String str = "drop "+(String)objectVector.elementAt(index);
					doparse(str);
					}
				}
			}
		});
	objectsList2.addFocusListener(new FocusAdapter()
		{
		public void focusGained(FocusEvent fe)
			{
			int index = objectsList2.getSelectedIndex();
			if(index == -1)
				{
				objectsList2.setSelectedIndex(objectsList2.getFirstVisibleIndex());
				}
			}

		});
	objectsScrollPane = new JScrollPane(objectsList2);
	objectsSplitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,objectsList1,objectsScrollPane);
	objectsSplitPane.setOneTouchExpandable(false);
	objectsSplitPane.setDividerSize(5);
	objectsSplitPane.setContinuousLayout(true);
	objects.addInternalFrameListener(new InternalFrameAdapter()
		{
		public void internalFrameActivated(InternalFrameEvent ife)
			{
			objectsList1.requestFocus();
			obj.setSelected(true);
			}
		});
	objects.getContentPane().add(objectsSplitPane,BorderLayout.CENTER);
	objectsSplitPane.setDividerLocation(100);
	jDesktopPanel.add(objects);
	objects.setBounds(registryArray[14],registryArray[15],registryArray[16],registryArray[17]);

//Einstellungen des Hauptfensters

	addWindowListener(this);
	setSize(registryArray[0],registryArray[1]);
	setTitle("BayesX");
	setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
	setIconImage(Toolkit.getDefaultToolkit().getImage("bayesicon.gif"));

// Random Number Generation

	Random seed = new Random();
	long[] seedArray = new long[6];
	for(int i=0; i<3; i++)
		{
		long help=(Long.valueOf("4294967087")).longValue();
		while(help>=(Long.valueOf("4294967087")).longValue() || help<=0)
			{
			help = seed.nextLong()/(Long.valueOf("10000000000")).longValue();
			}
		seedArray[i] = help;
		}
	for(int i=3; i<6; i++)		{
		long help=(Long.valueOf("4294944443")).longValue();
		while(help>=(Long.valueOf("4294944443")).longValue() || help<=0)
			{
			help = seed.nextLong()/(Long.valueOf("10000000000")).longValue();
			}
		seedArray[i] = help;
		}

	rStream = new RandMrg();
	rStream.setSeed(seedArray);
	}


//Die Listener zur Event-Bearbeitung:

//Listener fuer die Menueleiste und die Buttons

//ActionListener fuer File, die Buttons und die ComboBox

public void actionPerformed(ActionEvent ae)
	{
	String source = ae.getActionCommand();
	if(source.equals("Clear Output Window"))
		{
		fileCommand(2);
		}
	else if(source.equals("Open Output File"))
		{
		fileCommand(3);
		}
	else if(source.equals("Save Output"))
		{
		saveOutput();
		}
	else if(source.equals("Save Output as..."))
		{
		int returnVal = fileChooser.showSaveDialog(this);
		if(returnVal == JFileChooser.APPROVE_OPTION)
			{
			try
				{
				outputFile = fileChooser.getSelectedFile();
				if(getExtension(outputFile)!=null && (getExtension(outputFile).equals("rtf")||getExtension(outputFile).equals("txt")))
					{
					}
				else
					{
					if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
						{
						outputFile = new File(outputFile.getPath()+".rtf");
						}
					else if(fileChooser.getFileFilter().getDescription().equals("Text File (*.txt)"))
						{
						outputFile = new File(outputFile.getPath()+".txt");
						}
					else
						{
						}
					}
				if(outputFile.exists())
					{
					while(outputFile.exists())
						{
						Object[] options = {"Yes","No","Cancel"};
						int returnValue1 = JOptionPane.showOptionDialog(this, "File already exists!\nReplace existing file?",
							"Existing File", JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE,
							null, options, null);
						if(returnValue1==JOptionPane.YES_OPTION)
							{
							FileOutputStream outfile = new FileOutputStream(outputFile);
							if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
								{
								outputEditor.write(outfile,outputDocument,0,outputDocument.getLength());
								}
							else
								{
								printEditor.write(outfile,outputDocument,0,outputDocument.getLength());
								}
							outfile.close();
							hasBeenSaved = true;
							isSaved = true;
							save.setEnabled(false);
							return;
							}
						else
							{
							int returnValue2 = fileChooser.showSaveDialog(this);
							if(returnValue2 == JFileChooser.APPROVE_OPTION)
								{
								outputFile = fileChooser.getSelectedFile();
								if(getExtension(outputFile)!=null && (getExtension(outputFile).equals("rtf")||getExtension(outputFile).equals("txt")))
									{
									}
								else
									{
									if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
										{
										outputFile = new File(outputFile.getPath()+".rtf");
										}
									else if(fileChooser.getFileFilter().getDescription().equals("Text File (*.txt)"))
										{
										outputFile = new File(outputFile.getPath()+".txt");
										}
									else
										{
										}
									}
								}
							else
								{
								break;
								}
							}
						}
					FileOutputStream outfile = new FileOutputStream(outputFile);
					if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
						{
						outputEditor.write(outfile,outputDocument,0,outputDocument.getLength());
						}
					else
						{
						printEditor.write(outfile,outputDocument,0,outputDocument.getLength());
						}
					outfile.close();
					hasBeenSaved = true;
					isSaved = true;
					save.setEnabled(false);
					}
				else
					{
					FileOutputStream outfile = new FileOutputStream(outputFile);
					if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
						{
						outputEditor.write(outfile,outputDocument,0,outputDocument.getLength());
						}
					else
						{
						printEditor.write(outfile,outputDocument,0,outputDocument.getLength());
						}
					outfile.close();
					hasBeenSaved = true;
					isSaved = true;
					save.setEnabled(false);
					}
				}
			catch(IOException ioe)
				{
				System.err.println(ioe.getMessage());
				}
			catch(BadLocationException ble)
				{
				System.err.println(ble.getMessage());
				}
			}
		else
			{
			}
		}
	else if(source.equals("Print Output"))
		{
		try
			{
			FileOutputStream outfile = new FileOutputStream(printFile);
			printEditor.write(outfile,outputDocument,0,outputDocument.getLength());
			outfile.close();
			PrinterJob job = PrinterJob.getPrinterJob();
			PageFormat format = job.pageDialog(job.defaultPage());
			job.setPageable(new PageableText(printFile, format));
			if(job.printDialog())
				{
				job.print();
				}
			}
		catch(PrinterException pe)
			{
			System.err.println(pe.getMessage());
			}
		catch(IOException ioe)
			{
			System.err.println(ioe.getMessage());
			}
		catch(BadLocationException ble)
			{
			System.err.println(ble.getMessage());
			}
		}
	else if(source.equals("Exit"))
		{
		if(processRunning)
			{
//			parseThread.suspend();
			parseThread.setPriority(Thread.MIN_PRIORITY);
			setPause(true);
			Object[] options = {"Yes","No"};
			int returnVal = JOptionPane.showOptionDialog(this, "Are you sure?\nDo you want to interrupt the current process?",
				"Break", JOptionPane.YES_OPTION, JOptionPane.QUESTION_MESSAGE,
 				null, options, null);
			if(returnVal == JOptionPane.YES_OPTION)
				{
				parseThread.interrupt();
				setProcessrunning(false);
				finished();
//				setPause(false);
				pause = false;
				parseThread.setPriority(Thread.NORM_PRIORITY);
				pauseButton.setText("PAUSE");
				pauseButton.setToolTipText("Interrupts the current process");
				fileCommand(1);
				}
			else
				{
//				parseThread.resume();
				setPause(false);
				parseThread.setPriority(threadPriority);
				}
			}
		else
			{
			fileCommand(1);
			}
		}
	else if(source.equals("Default Windowing"))
		{
		this.setSize(750,500);
		command.setBounds(5,5,450,105);
		output.setBounds(5,115,450,285);
		review.setBounds(460,5,275,200);
		objects.setBounds(460,210,275,190);
		}
	else if(source.equals("Fontsize Small"))
		{
		scriptsize=10;
		fontFactor=0.84;
		StyleConstants.setFontSize(outputAttributes,scriptsize);
		commandArea.setFont(new Font("Monospaced",Font.PLAIN,10));
		reviewList.setFont(new Font("Sansserif",Font.PLAIN,10));
		objectsList1.setFont(new Font("Sansserif",Font.PLAIN,10));
		objectsList2.setFont(new Font("Sansserif",Font.PLAIN,10));
		}
	else if(source.equals("Fontsize Normal"))
		{
		scriptsize=12;
		fontFactor=1.0;
		StyleConstants.setFontSize(outputAttributes,scriptsize);
		commandArea.setFont(new Font("Monospaced",Font.PLAIN,12));
		reviewList.setFont(new Font("Sansserif",Font.PLAIN,12));
		objectsList1.setFont(new Font("Sansserif",Font.PLAIN,12));
		objectsList2.setFont(new Font("Sansserif",Font.PLAIN,12));
		}
 	else if(source.equals("Fontsize Large"))
		{
		scriptsize=16;
		fontFactor=1.4;
		StyleConstants.setFontSize(outputAttributes,scriptsize);
		commandArea.setFont(new Font("Monospaced",Font.PLAIN,16));
		reviewList.setFont(new Font("Sansserif",Font.PLAIN,16));
		objectsList1.setFont(new Font("Sansserif",Font.PLAIN,16));
		objectsList2.setFont(new Font("Sansserif",Font.PLAIN,16));
		}
         else if(source.equals("Reference Manual"))
                {
		Runtime r = Runtime.getRuntime();
		try
			{
			String[] comms = new String[2];
			comms[0]="ShellExec";
			comms[1]="doc\\reference_manual.pdf";
			Process p = r.exec(comms);
			}
		catch(IOException ioe)
			{
			System.err.println(ioe.getMessage());
			}
                }
         else if(source.equals("Methodology Manual"))
                {
		Runtime r = Runtime.getRuntime();
		try
			{
			String[] comms = new String[2];
			comms[0]="ShellExec";
			comms[1]="doc\\methodology_manual.pdf";
			Process p = r.exec(comms);
			}
		catch(IOException ioe)
			{
			System.err.println(ioe.getMessage());
			}
                }
         else if(source.equals("Tutorials Manual"))
                {
		Runtime r = Runtime.getRuntime();
		try
			{
			String[] comms = new String[2];
			comms[0]="ShellExec";
			comms[1]="doc\\tutorials_manual.pdf";
			Process p = r.exec(comms);
			}
		catch(IOException ioe)
			{
			System.err.println(ioe.getMessage());
			}
                }
        else if(source.equals("About BayesX"))
                {
		JOptionPane.showMessageDialog(this,"BayesX\n\nSoftware for Bayesian Inference in Structured Additive Regression Models\n"+
			"Version 2.1 (07.05.2012)\n\n"+
			"developed by\n"+
			"  Christiane Belitz\n"+
			"  Andreas Brezger\n"+
			"  Thomas Kneib\n"+
			"  Stefan Lang\n"+
			"  Nikolaus Umlauf\n\n"+
			"with contributions by\n"+
			"  Daniel Adler\n"+
			"  Eva-Maria Fronk\n"+
			"  Felix Heinzl\n"+
			"  Andrea Hennerfeind\n"+
			"  Manuela Hummel\n"+
			"  Alexander Jerak\n"+
			"  Susanne Konrath\n"+
			"  Petra Kragler\n"+
			"  Cornelia Oberhauser\n"+
			"  Leyre Estibaliz Osuna Echavarria\n"+
			"  Daniel Sabanes Bove\n"+
			"  Achim Zeileis\n\n"+
			"supported by\n"+
			"  Ludwig Fahrmeir (mentally)\n"+
			"  Leo Held (mentally)\n"+
			"  German Research Foundation (financially)","About BayesX"
			,JOptionPane.INFORMATION_MESSAGE,new ImageIcon(Toolkit.getDefaultToolkit().getImage("bayesicon.gif")));
		}
	else if(source.equals("BREAK"))
		{
		if(processRunning)
			{
//			parseThread.suspend();
			setPause(true);
			parseThread.setPriority(Thread.MIN_PRIORITY);
			Object[] options = {"Yes","No"};
			int returnVal = JOptionPane.showOptionDialog(this, "Are you sure?\nDo you want to interrupt the current process?",
				"Break", JOptionPane.YES_OPTION, JOptionPane.QUESTION_MESSAGE,
 				null, options, null);
			if(returnVal == JOptionPane.YES_OPTION)
				{
				setProcessrunning(false);
                                setStop(true);
				setPause(false);
				parseThread.setPriority(Thread.NORM_PRIORITY);
				pause = false;
				pauseButton.setText("PAUSE");
				pauseButton.setToolTipText("Interrupts the current process");
				Out("\nUSER BREAK\n\n");
				}
			else
				{
				if(!pause)
					{
//					parseThread.resume();
					setPause(false);
					parseThread.setPriority(threadPriority);
					}
				}
			}
		else
			{
			}
		}
	else if(source.equals("PAUSE"))
		{
		setPause(true);
		pause = true;
		parseThread.setPriority(Thread.MIN_PRIORITY);
		pauseButton.setText("CONTINUE");
		if(processRunning)
			{
			Out("\nPROGRAM PAUSED\nClick CONTINUE to proceed\n\n");
			}
		pauseButton.setToolTipText("Continues the current process");

		}
	else if(source.equals("CONTINUE"))
		{
		setPause(false);
		pause = false;
		parseThread.setPriority(Thread.NORM_PRIORITY);
		pauseButton.setText("PAUSE");
		if(processRunning)
			{
			Out("\nCONTINUED\n\n");
			}
		pauseButton.setToolTipText("Interrupts the current process");
		}
	else if(source.equals("SUPPRESS OUTPUT"))
		{
		setSuppressoutput(true);
		outputButton.setText("SHOW OUTPUT");
		outputButton.setToolTipText("Continues the documentation of output");
		}
	else if(source.equals("SHOW OUTPUT"))
		{
		setSuppressoutput(false);
		outputButton.setText("SUPPRESS OUTPUT");
		outputButton.setToolTipText("Suppresses the documentation of output");
		}
	else
		{
		JComboBox cb = (JComboBox)ae.getSource();
		source = (String)cb.getSelectedItem();
		if(source.equals(" VERY LOW"))
			threadPriority = 1;
		else if(source.equals(" LOW"))
			threadPriority = 3;
		else if(source.equals(" NORMAL"))
			threadPriority = 5;
		else if(source.equals(" HIGH"))
			threadPriority = 7;
		else if(source.equals(" VERY HIGH"))
			threadPriority = 9;
		if(processRunning && !pause)
			parseThread.setPriority(threadPriority);
		}
	}

//ItemListener fuer Window und Edit

public void itemStateChanged(ItemEvent ie)
	{
	Object source = ie.getItemSelectable();
	if(source==comm && comm.isSelected())
		{
		try
			{
			if(command.isIcon())
				{
				command.setIcon(false);
				}
			command.moveToFront();
			if(!command.isSelected())
				{
				command.setSelected(true);
				}
			}
		catch(PropertyVetoException pve)
			{
			System.err.println(pve.getMessage());
			}
		}
	else if(source==out && out.isSelected())
		{
		try
			{
			if(output.isIcon())
				{
				output.setIcon(false);
				}
			output.moveToFront();
			if(!output.isSelected())
				{
				output.setSelected(true);
				}
			}
		catch(PropertyVetoException pve)
			{
			System.err.println(pve.getMessage());
			}
		}
	else if(source==rev && rev.isSelected())
		{
		try
			{
			if(review.isIcon())
				{
				review.setIcon(false);
				}
			review.moveToFront();
			if(!review.isSelected())
				{
				review.setSelected(true);
				}
			}
		catch(PropertyVetoException pve)
			{
			System.err.println(pve.getMessage());
			}
		if(reviewList.getSelectedIndex()==-1 && reviewVector.size()>0)
			{
			reviewList.setSelectedIndex(reviewList.getFirstVisibleIndex());
			}
		}
	else if(source==obj && obj.isSelected())
		{
		try
			{
			if(objects.isIcon())
				{
				objects.setIcon(false);
				}
			objects.moveToFront();
			if(!objects.isSelected())
				{
				objects.setSelected(true);
				}
			}
		catch(PropertyVetoException pve)
			{
			System.err.println(pve.getMessage());
			}
		if(objectsList1.getSelectedIndex()==-1)
			{
			objectsList1.setSelectedIndex(objectsList1.getFirstVisibleIndex());
			}
		}
	}


//WindowListener

public void windowClosed(WindowEvent we)
	{
	}

public void windowOpened(WindowEvent we)
	{
	if(we.getSource() != mapFrame)
		{
		try
			{
			command.setSelected(true);
			StyleConstants.setFontSize(outputAttributes,scriptsize);
			commandArea.setFont(new Font("Monospaced",Font.PLAIN,scriptsize));
			reviewList.setFont(new Font("Sansserif",Font.PLAIN,scriptsize));
			objectsList1.setFont(new Font("Sansserif",Font.PLAIN,scriptsize));
			objectsList2.setFont(new Font("Sansserif",Font.PLAIN,scriptsize));
			if(scriptsize==10)
				{
				small.setSelected(true);
				}
			else if(scriptsize==12)
				{
				normal.setSelected(true);
				}
			else
				{
				large.setSelected(true);
				}
			}
		catch(PropertyVetoException pve)
			{
			System.err.println(pve.getMessage());
			}
		}
	}

public void windowIconified(WindowEvent we)
	{
	}

public void windowDeiconified(WindowEvent we)
	{
	}

public void windowActivated(WindowEvent we)
	{
	}

public void windowDeactivated(WindowEvent we)
	{
	}

public void windowClosing(WindowEvent we)
	{
	if(we.getSource()==mapFrame)
		{
		Object[] options = {"Yes","No","Cancel"};
		int returnVal = JOptionPane.showOptionDialog(mapFrame, "Save Graph?",
			"Save Graph", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE,
			null, options, null);
		if(returnVal == JOptionPane.YES_OPTION)
			{
			int returnVal1 = fileChooser2.showSaveDialog(mapFrame);
			if(returnVal1 == JFileChooser.APPROVE_OPTION)
				{
				try
					{
					File f = fileChooser2.getSelectedFile();
					if(getExtension(f)!=null && getExtension(f).equals("ps"))
						{
						}
					else
						{
						f = new File(f.getPath()+".ps");
						}
					if(f.exists())
						{
						while(f.exists())
							{
							Object[] options1 = {"Yes","No","Cancel"};
							int returnValue1 = JOptionPane.showOptionDialog(mapFrame, "File already exists!\nReplace existing file?",
								"Existing File", JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE,
								null, options1, null);
							if(returnValue1==JOptionPane.YES_OPTION)
								{
			                                        PrintWriter out = new PrintWriter(new FileWriter(f));
                                                                if(function==3)
                                                                        mapPanel.Saveplotnonp(out);
                                                                else if(function==2)
                                                                        mapPanel.Savedrawmap(out);
                                                                else if(function==1)
                                                                        mapPanel.SaveMap(out);
                                                                else if(function==4)
                                                                        mapPanel.Saveplotsample(out);
                                                                else if(function==5)
                                                                        mapPanel.Saveplotautocor(out);
                                                                else if(function==6)
                                                                        mapPanel.Saveplotsurf(out);
                                                                setEnabled(true);
                                                                mapFrame.dispose();
//								parseThread.resume();
								setPause(false);
								parseThread.setPriority(threadPriority);
                                                                return;
								}
							else
								{
								int returnValue2 = fileChooser2.showSaveDialog(mapFrame);
								if(returnValue2 == JFileChooser.APPROVE_OPTION)
									{
									f = fileChooser2.getSelectedFile();
									if(getExtension(f)!=null && getExtension(f).equals("ps"))
										{
										}
									else
										{
										f = new File(f.getPath()+".ps");
										}
									}
								else
									{
									break;
									}
								}
							}
	                                        PrintWriter out = new PrintWriter(new FileWriter(f));
                                                if(function==3)
                                                        mapPanel.Saveplotnonp(out);
                                                else if(function==2)
                                                        mapPanel.Savedrawmap(out);
                                                else if(function==1)
                                                        mapPanel.SaveMap(out);
                                                else if(function==4)
                                                        mapPanel.Saveplotsample(out);
                                                else if(function==5)
                                                        mapPanel.Saveplotautocor(out);
                                                else if(function==6)
                                                        mapPanel.Saveplotsurf(out);
                                                setEnabled(true);
						mapFrame.dispose();
//						parseThread.resume();
						setPause(false);
						parseThread.setPriority(threadPriority);
						}
					else
						{
                        	                PrintWriter out = new PrintWriter(new FileWriter(f));
                                                if(function==3)
                                                        mapPanel.Saveplotnonp(out);
                                                else if(function==2)
                                                        mapPanel.Savedrawmap(out);
                                                else if(function==1)
                                                        mapPanel.SaveMap(out);
                                                else if(function==4)
                                                        mapPanel.Saveplotsample(out);
                                                else if(function==5)
                                                        mapPanel.Saveplotautocor(out);
                                                else if(function==6)
                                                        mapPanel.Saveplotsurf(out);
                                                setEnabled(true);
						mapFrame.dispose();
//						parseThread.resume();
						setPause(false);
						parseThread.setPriority(threadPriority);
						}
					}
				catch(IOException ioe)
					{
					System.err.println(ioe.getMessage());
					}
				}
			}
		else if(returnVal==JOptionPane.NO_OPTION)
			{
 			setEnabled(true);
			mapFrame.dispose();
//			parseThread.resume();
			setPause(false);
			parseThread.setPriority(threadPriority);
			}
		else
			{
			}
		}
	else
		{
		if(processRunning)
			{
//			parseThread.suspend();
			setPause(true);
			parseThread.setPriority(Thread.MIN_PRIORITY);
			Object[] options = {"Yes","No"};
			int returnVal = JOptionPane.showOptionDialog(this, "Process running?\nDo you want to interrupt the current process?",
				"Break", JOptionPane.YES_OPTION, JOptionPane.QUESTION_MESSAGE,
					null, options, null);
			if(returnVal == JOptionPane.YES_OPTION)
				{
				parseThread.interrupt();
				setProcessrunning(false);
				finished();
//				setPause(false);
				pause = false;
				parseThread.setPriority(Thread.NORM_PRIORITY);
				pauseButton.setText("PAUSE");
				pauseButton.setToolTipText("Interrupts the current process");
				fileCommand(1);
				}
			else
				{
//				parseThread.resume();
				setPause(false);
				parseThread.setPriority(threadPriority);
				}
			}
		else
			{
			fileCommand(1);
			}
		}
	}


//ListSelectionListener fuer Review und Object

public void valueChanged(ListSelectionEvent lse)
	{
	if(lse.getSource() == objectsList1)
		{
		int index = objectsList1.getSelectedIndex();
		if(index>-1)
			{
			objectVector.clear();
			String selected = (String)objectsVector.elementAt(index);
			setObjectList(objectVector,selected);
			objectsList2.setListData(objectVector);
			}
		}
	}


//DocumentListener fuer das Output-Fenster

public void insertUpdate(DocumentEvent e)
	{
	isSaved = false;
	save.setEnabled(true);
	}

public void removeUpdate(DocumentEvent e)
	{
	isSaved = false;
	save.setEnabled(true);
	}

public void changedUpdate(DocumentEvent e)
	{
	isSaved = false;
	save.setEnabled(true);
	}



//AdjustmentListener fuer die Tabellen in JavaShowData

public void adjustmentValueChanged(AdjustmentEvent ae)
	{
	if(ae.getAdjustable()==horizontalScrollPane.getHorizontalScrollBar())
		{
		int c = horizontalScrollPane.getHorizontalScrollBar().getValue()/75;
		int r = verticalScrollBar.getValue();
		if(cols+c>varnames.size())
			{
			c=varnames.size()-cols;
			}
		if(rows+r>getRows())
			{
			r = getRows()-rows;
			}
		for(int i=0; i<cols; i++)
			{
			for(int j=0; j<rows; j++)
				{
				objectTable.setValueAt(getValue(j+r,i+c),j,i+c);
				}
			}
		}
	if(ae.getAdjustable()==verticalScrollBar)
		{
		int r = ae.getValue();
		int c = horizontalScrollPane.getHorizontalScrollBar().getValue()/75;
		if(cols+c>varnames.size())
			{
			c=varnames.size()-cols;
			}
		if(rows+r>getRows())
			{
			r = getRows()-rows;
			}
		for(int j=0; j<rows; j++)
			{
			for(int i=0; i<cols; i++)
				{
				objectTable.setValueAt(getValue(j+r,i+c),j,i+c);
				}
			rowHeaderTable.setValueAt(new Integer(j+r+1),j+1,0);
			}
		}
	}

//Hilfsfunktionen zur Event-Bearbeitung

//liefert die Endung einer Datei zurueck

public String getExtension(File f)
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

//Zeigt ein Dataset in einem separaten Fenster an


public void JavaShowData()
	{
	setEnabled(false);

	varnames = new Vector();
	setVarnames(varnames);
	rows = Math.min(19,getRows());
	cols = Math.min(7,varnames.size());
	data = new Object[rows][varnames.size()];
	vars = new Object[varnames.size()];

	for(int i=0; i<cols; i++)
		{
		for(int j=0; j<rows; j++)
			{
			data[j][i]=getValue(j,i);
			}
		}
	for(int i=0; i<varnames.size();i++)
		{
		vars[i] = varnames.elementAt(i);
		}
	tableModel1 = new DefaultTableModel(data,vars);
	objectTable = new JTable(tableModel1);
	objectTable.setFont(new Font("Monospaced",Font.PLAIN,scriptsize));
	objectTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
	objectTable.setPreferredScrollableViewportSize(new Dimension(cols*75,rows*16));
	objectTable.setEnabled(false);
	horizontalScrollPane = new JScrollPane(objectTable, JScrollPane.VERTICAL_SCROLLBAR_NEVER, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
	horizontalScrollPane.getHorizontalScrollBar().setUnitIncrement(75);
	horizontalScrollPane.getHorizontalScrollBar().setBlockIncrement(75*(cols-1));
	horizontalScrollPane.getHorizontalScrollBar().addAdjustmentListener(this);

	verticalScrollBar = new JScrollBar();
	verticalScrollBar.setMaximum(getRows());
	verticalScrollBar.setBlockIncrement(rows-1);
	verticalScrollBar.setVisibleAmount(objectTable.getRowCount());
	verticalScrollBar.addAdjustmentListener(this);

	rowHeaderData1 = new Object[rows+1][1];
	rowHeaderData2 = new Object[1];
	for(int i=1; i<rows+1; i++)
		{
		rowHeaderData1[i][0] = new Integer(i);
		}
	rowHeaderData2[0] = "";
	tableModel2 = new DefaultTableModel(rowHeaderData1, rowHeaderData2);
	rowHeaderTable = new JTable(tableModel2);
	rowHeaderTable.setFont(new Font("Monospaced",Font.PLAIN,scriptsize));
	rowHeaderTable.setBackground(Color.getColor("lightgray"));
	rowHeaderTable.setEnabled(false);

	objectFrame = new JFrame("Object -Viewer");
	objectFrame.setIconImage(Toolkit.getDefaultToolkit().getImage("bayesicon.gif"));
	objectFrame.getContentPane().add(horizontalScrollPane,BorderLayout.CENTER);
	objectFrame.getContentPane().add(verticalScrollBar,BorderLayout.EAST);
	objectFrame.getContentPane().add(rowHeaderTable,BorderLayout.WEST);
 	objectFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
	objectFrame.addWindowListener(new WindowAdapter()
		{
		public void windowClosing(WindowEvent evt)
			{
			setEnabled(true);
			objectFrame.dispose();
//			parseThread.resume();
			setPause(false);
			parseThread.setPriority(threadPriority);
			}
		});
	objectFrame.addComponentListener(new ComponentAdapter()
		{
		public void componentResized(ComponentEvent ce)
			{
			int r = verticalScrollBar.getValue();
			int c = horizontalScrollPane.getHorizontalScrollBar().getValue()/75;
			cols = horizontalScrollPane.getHorizontalScrollBar().getVisibleAmount()/75;
			if(horizontalScrollPane.getHorizontalScrollBar().getVisibleAmount()%75>0)
				{
				cols++;
				}
			cols = Math.min(cols,varnames.size());
			horizontalScrollPane.getHorizontalScrollBar().setBlockIncrement(75*(cols-1));
			if(cols+c>varnames.size())
				{
				c=varnames.size()-cols;
				}
			int newRows = objectFrame.getHeight()/16-3;
			newRows = Math.min(newRows,getRows());
			if(newRows+r>getRows())
				{
				r = getRows()-newRows;
				}
			if(newRows!=rows)
				{
				rows = Math.min(newRows,getRows());
				data = new Object[rows][varnames.size()];
				rowHeaderData1 = new Object[rows+1][1];
				tableModel1.setDataVector(data, vars);
				tableModel2.setDataVector(rowHeaderData1, rowHeaderData2);
				verticalScrollBar.setVisibleAmount(rows);
				verticalScrollBar.setBlockIncrement(rows-1);
				}
			for(int j=0; j<rows; j++)
				{
				for(int i=0; i<cols; i++)
					{
					objectTable.setValueAt(getValue(j+r,i+c),j,i+c);
					}
				rowHeaderTable.setValueAt(new Integer(j+r+1),j+1,0);
				}
			}
		});
	objectFrame.pack();
	objectFrame.show();
//	parseThread.suspend();
	setPause(true);
	parseThread.setPriority(Thread.MIN_PRIORITY);
	}

//Funktionen zum Zeichen von Karten

public void JavaDescribeMap(boolean opt)
	{

        int width;
        int height;

        function = 1;

        drawnames = opt;
	fontsize = 0;
        outfile = "";

        setEnabled(false);
	fileChooser2 = new JFileChooser(defaultDirectory);
	fileChooser2.addChoosableFileFilter(new MapFilter());
	mapFrame = new JFrame("Object-Viewer");
	double[] d = new double[4];
	getboundaries(d);
        mapPanel = new MapPanel(this);

//        mapPanel.setfunction((short)1);

	mapFrame.getContentPane().add(mapPanel,BorderLayout.CENTER);
	mapFrame.setIconImage(Toolkit.getDefaultToolkit().getImage("bayesicon.gif"));
 	mapFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
 	mapFrame.addWindowListener(this);

        width = 500;
        height = (int)(width*(d[3]-d[1])/(d[2]-d[0]));

        if(height>700)
                {
                width = (int)(width*700/height);
                height = 700;
                }

        mapFrame.setSize(width + 8,height + 27);


        if(outfile.equals(""))
                {
		mapFrame.show();
//		parseThread.suspend();
		setPause(true);
		parseThread.setPriority(Thread.MIN_PRIORITY);
		}
        else
                {
                try
                        {
                        PrintWriter out = new PrintWriter(new FileWriter(outfile));
                        mapPanel.SaveMap(out);
                        setEnabled(true);
                        mapFrame.dispose();
                        }
                catch(IOException ioe)
        	        {
	        	System.err.println(ioe.getMessage());
		        }
                }

	}

public void JavaShowMap(boolean opt, int jfontsize, String joutfile)
	{

        int width;
        int height;

        function = 1;

        drawnames = opt;
	fontsize = jfontsize;
        outfile = joutfile;

        setEnabled(false);
	fileChooser2 = new JFileChooser(defaultDirectory);
	fileChooser2.addChoosableFileFilter(new MapFilter());
	mapFrame = new JFrame("Object-Viewer");
	double[] d = new double[4];
	getboundaries(d);
        mapPanel = new MapPanel(this);

//        mapPanel.setfunction((short)1);

	mapFrame.getContentPane().add(mapPanel,BorderLayout.CENTER);
	mapFrame.setIconImage(Toolkit.getDefaultToolkit().getImage("bayesicon.gif"));
 	mapFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
 	mapFrame.addWindowListener(this);

        width = 500;
        height = (int)(width*(d[3]-d[1])/(d[2]-d[0]));

        if(height>700)
                {
                width = (int)(width*700/height);
                height = 700;
                }

        mapFrame.setSize(width + 8,height + 27);


        if(outfile.equals(""))
                {
		mapFrame.show();
//		parseThread.suspend();
		setPause(true);
		parseThread.setPriority(Thread.MIN_PRIORITY);
		}
        else
                {
                try
                        {
                        PrintWriter out = new PrintWriter(new FileWriter(outfile));
                        mapPanel.SaveMap(out);
                        setEnabled(true);
                        mapFrame.dispose();
                        }
                catch(IOException ioe)
        	        {
	        	System.err.println(ioe.getMessage());
		        }
                }

	}

public void Javadrawmap(boolean opt1, boolean opt2, boolean opt3, boolean opt4, boolean opt5,
                        double jlowerlimit, double jupperlimit, short jshades, boolean jpcat,
                        int jfontsize, String joutfile, String jtitle, double jtitlescale)
	{
        if(joutfile.equals(""))
                {
//                parseThread.suspend();
		setPause(true);
		parseThread.setPriority(Thread.MIN_PRIORITY);
		}

        int width;
        int height;

        function = 2;

        color = opt1;
        legend = opt2;
        swap = opt3;
        drawnames = opt4;
	hcl = opt5;
        lowerlimit = jlowerlimit;
        upperlimit = jupperlimit;
        shades = jshades;
	fontsize = jfontsize;
        outfile = joutfile;
        title = jtitle;
        titlescale = jtitlescale;

        if(jpcat)
            {
            shades = 3;
            lowerlimit = -1;
            upperlimit = 1;
            }

        setEnabled(false);
	double[] d = new double[4];
        getboundaries(d);
	mapFrame = new JFrame("Object-Viewer");
	mapPanel = new MapPanel(this);

//        mapPanel.setfunction((short)2);

        if(outfile.equals(""))
                {
	        fileChooser2 = new JFileChooser(defaultDirectory);
	        fileChooser2.addChoosableFileFilter(new MapFilter());
                mapFrame.getContentPane().add(mapPanel,BorderLayout.CENTER);
	        mapFrame.setIconImage(Toolkit.getDefaultToolkit().getImage("bayesicon.gif"));
 	        mapFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
 	        mapFrame.addWindowListener(this);
                }

        width = 500;
        height = (int)(width*(d[3]-d[1])/(d[2]-d[0]));

        int help = 0;
        if(legend)
                {
                height = height + 45;
                help = help + 45;
                }
        if(!title.equals(""))
                {
                height = height + 35;
                help = help + 35;
                }

        if(height>700)
                {
                width = (int)(width*(700-help)/height);
                height = 700;
                }

        mapFrame.setSize(width + 8, height + 27);

        if(outfile.equals(""))
                {
                mapFrame.show();
                }
        else
                {
                try
                        {
                        PrintWriter out = new PrintWriter(new FileWriter(outfile));
                        mapPanel.Savedrawmap(out);
                        setEnabled(true);
                        mapFrame.dispose();
                        }
                catch(IOException ioe)
        	        {
	        	System.err.println(ioe.getMessage());
		        }
                }

                if(nrNA > 0 && nrNA < getnrregions())
                    {
                    String str = "NOTE: " + String.valueOf(nrNA) + " missing value(s) plotted\n";
                    Out(str,false,false,(short)11,0,0,0);
                    }
                else if (nrNA >= getnrregions())
                    {
                    String str = "WARNING: only missing values plotted - map probably doesn't match data file\n";
                    Out(str,true,true,(short)11,255,0,0);
                    }

	}

public void Javaplotnonp(String joutfile, String jtitle, String jxlab, String jylab,
                         String jconnect, String jlinecolor, int jheight, int jwidth, double jxmax, double jxmin,
                         double jymax, double jymin, double jxstep, double jxstart, double jystep, double jystart,
			 int jyear, int jmonth, int jlinewidth, int jpointsize, int jfontsize, double jtitlescale)

        {
        if(joutfile.equals(""))
                {
//                parseThread.suspend();
		setPause(true);
		parseThread.setPriority(Thread.MIN_PRIORITY);
		}
        function = 3;

        outfile = joutfile;
        title = jtitle;
        xlab = jxlab;
        ylab = jylab;
        connect = jconnect;
	linecolor = jlinecolor;
        height = jheight;
        width = jwidth;
        xmax = jxmax;
        xmin = jxmin;
        ymax = jymax;
        ymin = jymin;
        xstep = jxstep;
        ystep = jystep;
        xstart = jxstart;
        ystart = jystart;
        year = jyear;
        month = jmonth;
	linewidth = jlinewidth;
	pointsize = jpointsize;
	fontsize = jfontsize;
	titlescale = jtitlescale;

        setEnabled(false);

        mapFrame = new JFrame("Object-Viewer");

        mapPanel = new MapPanel(this);
//        mapPanel.setfunction((short)3);

//        MapPanel mapPanel2 = new MapPanel(this);
//        mapPanel2.setfunction((short)3);

        if(outfile.equals(""))
                {
	        fileChooser2 = new JFileChooser(defaultDirectory);
	        fileChooser2.addChoosableFileFilter(new MapFilter());

                JTabbedPane mapTabbedPane = new JTabbedPane();
	        mapTabbedPane.addTab("Page 1",mapPanel);
//                mapTabbedPane.addTab("Page 2",mapPanel2);

                mapFrame.getContentPane().add(mapTabbedPane,BorderLayout.CENTER);
//	        mapFrame.getContentPane().add(mapPanel,BorderLayout.CENTER);

	        mapFrame.setIconImage(Toolkit.getDefaultToolkit().getImage("bayesicon.gif"));
 	        mapFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
 	        mapFrame.addWindowListener(this);
                mapFrame.setSize(596,(height>210)?(height+240):500);
                mapFrame.show();
                }
        else
                {
                try
                        {
                        PrintWriter out = new PrintWriter(new FileWriter(outfile));
                        mapPanel.Saveplotnonp(out);
                        setEnabled(true);
                        mapFrame.dispose();
                        }
                catch(IOException ioe)
        	        {
	        	System.err.println(ioe.getMessage());
		        }
                }

//        resetplotoptions();

        }

/*
public void Javaplot(String joutfile, String jtitle, String jxlab, String jylab,
                     String jconnect, int jheight, int jwidth, int jlinewidth, int jfontsize)
        {
        function = 3;

        outfile = joutfile;
        title = jtitle;
        xlab = jxlab;
        ylab = jylab;
        connect = jconnect;
        height = jheight;
        width = jwidth;
	linewidth = jlinewidth;
	fontsize = jfontsize;

        setEnabled(false);

        mapFrame = new JFrame("Object-Viewer");

        if(outfile.equals(""))
                {
	        fileChooser2 = new JFileChooser(defaultDirectory);
	        fileChooser2.addChoosableFileFilter(new MapFilter());

                int nrpages = 5;    // berechnen

                JTabbedPane mapTabbedPane = new JTabbedPane();
                MapPanel[] mapPanels = new MapPanel[nrpages];

                for(int i=0;i<nrpages;i++)
                    {
                    mapPanels[i] = new MapPanel(this);
//                    mapPanels[i].setfunction((short)3);
// hier alles speichern, was man fuer die ganze Seite braucht
//                  for(int j=0;j<nrplots;j++)
//                  mapPanels[i].data[j][][]
//                  mapPanels[i].minX[j]
//                  mapPanels[i].setparam[i]

                    mapTabbedPane.addTab("Page "+(i+1),mapPanels[i]);
                    }

                mapFrame.getContentPane().add(mapTabbedPane,BorderLayout.CENTER);
	        mapFrame.setIconImage(Toolkit.getDefaultToolkit().getImage("bayesicon.gif"));
 	        mapFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
 	        mapFrame.addWindowListener(this);
                mapFrame.setSize(596,(height>210)?(height+240):500);
                mapFrame.show();

                parseThread.suspend();
                }
        else
                {
                try
                        {
                        PrintWriter out = new PrintWriter(new FileWriter(outfile));
                        mapPanel.Saveplotnonp(out);
                        setEnabled(true);
                        mapFrame.dispose();
                        }
                catch(IOException ioe)
        	        {
	        	System.err.println(ioe.getMessage());
		        }
                }
        }
*/

public void Javaplotsample(String joutfile, String jconnect)
        {
        if(joutfile.equals(""))
                {
//                parseThread.suspend();
		setPause(true);
		parseThread.setPriority(Thread.MIN_PRIORITY);
		}
        function = 4;

        plotsperpage = 6;

        outfile = joutfile;
        connect = jconnect;

        setEnabled(false);
	mapFrame = new JFrame("Object-Viewer");

        mapPanel = new MapPanel(this);

        if(outfile.equals(""))
                {
	        fileChooser2 = new JFileChooser(defaultDirectory);
	        fileChooser2.addChoosableFileFilter(new MapFilter());

                int nrpages = (getDCols()-2)/plotsperpage+1;

        	JTabbedPane mapTabbedPane = new JTabbedPane();
                MapPanel[] mapPanels = new MapPanel[nrpages];

                for(int i=0;i<nrpages;i++)
                    {
                    mapPanels[i] = new MapPanel(this);
                    mapPanels[i].page = (short)(i);
                    mapTabbedPane.addTab("Page "+(i+1),mapPanels[i]);
                    }

                mapFrameScrollPane = new JScrollPane(mapTabbedPane);
//                mapFrameScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
                mapFrameScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);

//                mapFrame.getContentPane().add(mapTabbedPane,BorderLayout.CENTER);
                mapFrame.getContentPane().add(mapFrameScrollPane,BorderLayout.CENTER);
	        mapFrame.setIconImage(Toolkit.getDefaultToolkit().getImage("bayesicon.gif"));
 	        mapFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
 	        mapFrame.addWindowListener(this);
                mapFrame.setSize(596,740);
                mapFrame.show();
                }
        else
                {
                try
                        {
                        PrintWriter out = new PrintWriter(new FileWriter(outfile));
//                        mapPanel = new MapPanel(this);
                        mapPanel.Saveplotsample(out);
                        setEnabled(true);
                        mapFrame.dispose();
                        }
                catch(IOException ioe)
        	        {
	        	System.err.println(ioe.getMessage());
		        }
                }
        }

public void Javaplotautocor(String joutfile, String jconnect, boolean jmeanautocor)
        {
        if(joutfile.equals(""))
                {
//                parseThread.suspend();
		setPause(true);
		parseThread.setPriority(Thread.MIN_PRIORITY);
		}
        function = 5;

        outfile = joutfile;
        connect = jconnect;
        meanautocor = jmeanautocor;

        if(meanautocor==false)
            plotsperpage = 6;
        else
            plotsperpage = 3;

        setEnabled(false);
	mapFrame = new JFrame("Object-Viewer");

        mapPanel = new MapPanel(this);

        if(outfile.equals(""))
                {
	        fileChooser2 = new JFileChooser(defaultDirectory);
	        fileChooser2.addChoosableFileFilter(new MapFilter());

                int nrpages = (getDCols()-2)/plotsperpage+1;

        	JTabbedPane mapTabbedPane = new JTabbedPane();
                MapPanel[] mapPanels = new MapPanel[nrpages];

                for(int i=0;i<nrpages;i++)
                    {
                    mapPanels[i] = new MapPanel(this);
                    mapPanels[i].page = (short)(i);
                    mapTabbedPane.addTab("Page "+(i+1),mapPanels[i]);
                    }

                mapFrameScrollPane = new JScrollPane(mapTabbedPane);
                mapFrameScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);

//                mapFrame.getContentPane().add(mapTabbedPane,BorderLayout.CENTER);
                mapFrame.getContentPane().add(mapFrameScrollPane,BorderLayout.CENTER);
	        mapFrame.setIconImage(Toolkit.getDefaultToolkit().getImage("bayesicon.gif"));
 	        mapFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
 	        mapFrame.addWindowListener(this);
                mapFrame.setSize(596,740);
                mapFrame.show();
                }
        else
                {
                try
                        {
                        PrintWriter out = new PrintWriter(new FileWriter(outfile));
                        mapPanel.Saveplotautocor(out);
                        setEnabled(true);
                        mapFrame.dispose();
                        }
                catch(IOException ioe)
        	        {
	        	System.err.println(ioe.getMessage());
		        }
                }
        }


public void Javaplotsurf(String joutfile, String jtitle, String jxlab, String jylab, String jzlab,
			 double jxrot, double jyrot, double jzrot, String jlinecolor, int jheight, int jwidth,
 			 double jxlimtop, double jxlimbottom, double jylimtop, double jylimbottom, double jzlimtop, double jzlimbottom,
	                 double jxstep, double jxstart, double jystep, double jystart, double jzstep, double jzstart, int jgridsize,
			 int jlinewidth, int jpointsize, int jfontsize, double jtitlescale)

        {
        if(joutfile.equals(""))
                {
		setPause(true);
		parseThread.setPriority(Thread.MIN_PRIORITY);
		}
        function = 6;

        outfile = joutfile;
        title = jtitle;
        xlab = jxlab;
        ylab = jylab;
 	zlab = jzlab;
	xrot = jxrot;
	yrot = jyrot;
	zrot = jzrot;
	linecolor = jlinecolor;
        height = jheight;
        width = jwidth;
	xmax = jxlimtop;
	ymax = jylimtop;
	zmax = jzlimtop;
	xmin = jxlimbottom;
	ymin = jylimbottom;
	zmin = jzlimbottom;
	xstep = jxstep;
	ystep = jystep;
	zstep = jzstep;
	xstart = jxstart;
	ystart = jystart;
	zstart = jzstart;
	gridsize = jgridsize;
        linewidth = jlinewidth;
	pointsize = jpointsize;
	fontsize = jfontsize;
	titlescale = jtitlescale;

        setEnabled(false);

        mapFrame = new JFrame("Object-Viewer");

        mapPanel = new MapPanel(this);

        if(outfile.equals(""))
                {
	        fileChooser2 = new JFileChooser(defaultDirectory);
	        fileChooser2.addChoosableFileFilter(new MapFilter());

                mapFrame.getContentPane().add(mapPanel,BorderLayout.CENTER);
	        mapFrame.setIconImage(Toolkit.getDefaultToolkit().getImage("bayesicon.gif"));
 	        mapFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
 	        mapFrame.addWindowListener(this);
                mapFrame.setSize(596,(height>210)?(height+240):500);
                mapFrame.show();
                }
        else
                {
                try
                        {
                        PrintWriter out = new PrintWriter(new FileWriter(outfile));
                        mapPanel.Saveplotsurf(out);
                        setEnabled(true);
                        mapFrame.dispose();
                        }
                catch(IOException ioe)
        	        {
	        	System.err.println(ioe.getMessage());
		        }
                }

        }


/*
private void resetplotoptions()
        {
        xstep = 0.0;
        ystep = 0.0;
        }
*/

// fuer maps
native void getline(double[] d, int i, int j, int k);
native void getboundaries(double[] d);
native int getnrregions();
native int getnrpoly(int i);
native int getnrlines(int i,int j);
native boolean isin(int i);
native String getregionname(int i);
native void getcentroid(double[] centroid, int i);

native double getname(int i);

//Zwei Funktionen zum Steuern des Verhaltens beim Speichern bzw. Schliessen der Anwendung

public void fileAction(int i)
	{
	if(i==1)
		{
		try
			{
			if(registryFile.canWrite())
				{
				registryString = getWidth()+" "+getHeight()+" "+output.getX()+" "+output.getY()+" "+output.getWidth()+" "+output.getHeight()+" "+command.getX()+" "+command.getY()+
					" "+command.getWidth()+" "+command.getHeight()+" "+review.getX()+" "+review.getY()+" "+review.getWidth()+" "+review.getHeight()+" "+objects.getX()+
					" "+objects.getY()+" "+objects.getWidth()+" "+objects.getHeight()+" "+scriptsize+"\nBayesX System-File. DO NOT EDIT!";
				PrintWriter out = new PrintWriter(new FileWriter(registryFile));
				if(out!=null)
					{
					out.println(registryString);
					out.close();
					}
				}
			else
				{
//				System.out.println("write chack failed 4");
				}
			System.exit(0);
			}
		catch(IOException ioe)
			{
			System.err.println(ioe.getMessage());
			System.exit(0);
			}
		}
	else if(i==2)
		{
		outputPane.setText("");
		isSaved = true;
		hasBeenSaved = false;
		save.setEnabled(false);
		}
	else if(i==3)
		{
		try
			{
			int returnVal = fileChooser.showOpenDialog(this);
			if(returnVal == JFileChooser.APPROVE_OPTION)
				{
				outputPane.setText("");
				outputFile = fileChooser.getSelectedFile();
				FileInputStream infile = new FileInputStream(outputFile);

				if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
					{
					outputEditor.read(infile,outputDocument,0);
					}
				else
					{
					printEditor.read(infile,outputDocument,0);
					}
				infile.close();
				isSaved = true;
				hasBeenSaved = true;
				save.setEnabled(false);
				}
			}
		catch(IOException ioe)
			{
			System.err.println(ioe.getMessage());
			}
		catch(BadLocationException ble)
			{
			System.err.println(ble.getMessage());
			}
		}
	}


public void fileCommand(int i)
	{
	try
		{
		if(isSaved)
			{
			fileAction(i);
			}
		else
			{
			Object[] options = {"Yes","No","Cancel"};
			int returnVal = JOptionPane.showOptionDialog(this, "Output has changed!\nSave changes?",
				"Save Output", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE,
					null, options, null);
			if(returnVal == JOptionPane.YES_OPTION)
				{
				if(hasBeenSaved)
					{
					FileOutputStream outfile = new FileOutputStream(outputFile);
					if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
						{
						outputEditor.write(outfile,outputDocument,0,outputDocument.getLength());
						}
					else
						{
						printEditor.write(outfile,outputDocument,0,outputDocument.getLength());
						}
					outfile.close();
					fileAction(i);
					}
				else
					{
					int returnValue = fileChooser.showSaveDialog(this);
					if(returnValue == JFileChooser.APPROVE_OPTION)
						{
						outputFile = fileChooser.getSelectedFile();
						if(getExtension(outputFile)!=null && (getExtension(outputFile).equals("rtf")||getExtension(outputFile).equals("txt")))
							{
							}
						else
							{
							if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
								{
								outputFile = new File(outputFile.getPath()+".rtf");
								}
							else if(fileChooser.getFileFilter().getDescription().equals("Text File (*.txt)"))
								{
								outputFile = new File(outputFile.getPath()+".txt");
								}
							else
								{
								}
							}
						if(outputFile.exists())
							{
							while(outputFile.exists())
								{
								int returnValue1 = JOptionPane.showOptionDialog(this, "File already exists!\nReplace existing file?",
									"Existing File", JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE,
										null, options, null);
								if(returnValue1==JOptionPane.YES_OPTION)
									{
									FileOutputStream outfile = new FileOutputStream(outputFile);
									if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
										{
										outputEditor.write(outfile,outputDocument,0,outputDocument.getLength());
										}
									else
										{
										printEditor.write(outfile,outputDocument,0,outputDocument.getLength());
										}
									outfile.close();
									fileAction(i);
									return;
									}
								else
									{
									int returnValue2 = fileChooser.showSaveDialog(this);
									if(returnValue2 == JFileChooser.APPROVE_OPTION)
										{
										outputFile = fileChooser.getSelectedFile();
										if(getExtension(outputFile)!=null && (getExtension(outputFile).equals("rtf")||getExtension(outputFile).equals("txt")))
											{
											}
										else
											{
											if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
												{
												outputFile = new File(outputFile.getPath()+".rtf");
												}
											else if(fileChooser.getFileFilter().getDescription().equals("Text File (*.txt)"))
												{
												outputFile = new File(outputFile.getPath()+".txt");
												}
											else
												{
												}
											}
										}
									else
										{
										break;
										}
									}
								}
							FileOutputStream outfile = new FileOutputStream(outputFile);
							if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
								{
								outputEditor.write(outfile,outputDocument,0,outputDocument.getLength());
								}
							else
								{
								printEditor.write(outfile,outputDocument,0,outputDocument.getLength());
								}
							outfile.close();
							fileAction(i);
							}
						else
							{
							FileOutputStream outfile = new FileOutputStream(outputFile);
							if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
								{
								outputEditor.write(outfile,outputDocument,0,outputDocument.getLength());
								}
							else
								{
								printEditor.write(outfile,outputDocument,0,outputDocument.getLength());
								}
							outfile.close();
							fileAction(i);
							}
						}
					else
						{
						}
					}
				}
			else if(returnVal == JOptionPane.NO_OPTION)
				{
				fileAction(i);
				}
			else
				{
				}
			}
		}
	catch(IOException ioe)
		{
		System.err.println(ioe.getMessage());
		}
	catch(BadLocationException ble)
		{
		System.err.println(ble.getMessage());
		}
	}


//Funktion zum Speichern des Outputs

public void saveOutput()
	{
	try
		{
		if(hasBeenSaved)
			{
			FileOutputStream outfile = new FileOutputStream(outputFile);
			if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
				{
				outputEditor.write(outfile,outputDocument,0,outputDocument.getLength());
				}
			else
				{
				printEditor.write(outfile,outputDocument,0,outputDocument.getLength());
				}
			outfile.close();
			isSaved = true;
			save.setEnabled(false);
			}
		else
			{
			int returnVal = fileChooser.showSaveDialog(this);
			if(returnVal == JFileChooser.APPROVE_OPTION)
				{
				outputFile = fileChooser.getSelectedFile();
				if(getExtension(outputFile)!=null && (getExtension(outputFile).equals("rtf")||getExtension(outputFile).equals("txt")))
					{
					}
				else
					{
					if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
						{
						outputFile = new File(outputFile.getPath()+".rtf");
						}
					else if(fileChooser.getFileFilter().getDescription().equals("Text File (*.txt)"))
						{
						outputFile = new File(outputFile.getPath()+".txt");
						}
					else
						{
						}
					}
				if(outputFile.exists())
					{
					while(outputFile.exists())
						{
						Object[] options = {"Yes","No","Cancel"};
						int returnValue1 = JOptionPane.showOptionDialog(this, "File already exists!\nReplace existing file?",
							"Existing File", JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE,
							null, options, null);
						if(returnValue1==JOptionPane.YES_OPTION)
							{
							FileOutputStream outfile = new FileOutputStream(outputFile);
							if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
								{
								outputEditor.write(outfile,outputDocument,0,outputDocument.getLength());
								}
							else
								{
								printEditor.write(outfile,outputDocument,0,outputDocument.getLength());
								}
							outfile.close();
							hasBeenSaved = true;
							isSaved = true;
							save.setEnabled(false);
							return;
							}
						else
							{
							int returnValue2 = fileChooser.showSaveDialog(this);
							if(returnValue2 == JFileChooser.APPROVE_OPTION)
								{
								outputFile = fileChooser.getSelectedFile();
								if(getExtension(outputFile)!=null && (getExtension(outputFile).equals("rtf")||getExtension(outputFile).equals("txt")))
									{
									}
								else
									{
									if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
										{
										outputFile = new File(outputFile.getPath()+".rtf");
										}
									else if(fileChooser.getFileFilter().getDescription().equals("Text File (*.txt)"))
										{
										outputFile = new File(outputFile.getPath()+".txt");
										}
									else
										{
										}
									}
								}
							else
								{
								break;
								}
							}
						}
					FileOutputStream outfile = new FileOutputStream(outputFile);
					if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
						{
						outputEditor.write(outfile,outputDocument,0,outputDocument.getLength());
						}
					else
						{
						printEditor.write(outfile,outputDocument,0,outputDocument.getLength());
						}
					outfile.close();
					hasBeenSaved = true;
					isSaved = true;
					save.setEnabled(false);
					}
				else
					{
					FileOutputStream outfile = new FileOutputStream(outputFile);
					if(fileChooser.getFileFilter().getDescription().equals("Rich Text Format (*.rtf)"))
						{
						outputEditor.write(outfile,outputDocument,0,outputDocument.getLength());
						}
					else
						{
						printEditor.write(outfile,outputDocument,0,outputDocument.getLength());
						}
					outfile.close();
					hasBeenSaved = true;
					isSaved = true;
					save.setEnabled(false);
					}
				}
			else
				{
				}
			}
		}
	catch(IOException ioe)
		{
		System.err.println(ioe.getMessage());
		}
	catch(BadLocationException ble)
		{
		System.err.println(ble.getMessage());
		}
	}


private void checkDefaultpath(String str)
	{
	if(!(new File(str)).exists())
		{
		Out("ERROR: Path "+str+" does not exist.",true,true,(short)11,255,0,0);
		}
	else if(!(new File(str)).canWrite())
		{
		Out("WARNING: No permission to write to default directories.",true,true,(short)11,255,0,0);
		Out("         Specify a new default directory using the defaultpath command.\n"+
		    "         Type for example: defaultpath=c:\\temp",false,false,(short)11,0,0,0);
		}
	else
		{
		if(!(new File(str,"output")).exists())
			{
			(new File(str,"output")).mkdir();
			}
		if(!(new File(str,"temp")).exists())
			{
			(new File(str,"temp")).mkdir();
			}
		}
	}

//private void JavaOutput(String str)
//private void JavaOutput(String str, boolean thick, boolean italic, short size)
//public void JavaOutput(String str, boolean thick, boolean italic, short size, int c1, int c2, int c3)
private void JavaOutput(String str, boolean thick, boolean italic, short size, int c1, int c2, int c3)
        {
	if(consoleInput)
		{
		System.out.print(str);
		}
	else
		{
		try
			{
			final String stri = str;
			final boolean b1 = thick;
			final boolean b2 = italic;
			final int s = (int)(size*fontFactor);
			final int i1 = c1;
			final int i2 = c2;
			final int i3 = c3;
			SwingUtilities.invokeAndWait(new Runnable()
				{
				public void run()
					{
					try
						{
						if(i1!=0 || i2!=0 || i3!=0)
							{
							StyleConstants.setForeground(outputAttributes, new Color(i1, i2, i3));
							}
						if(b1==true)
							{
							StyleConstants.setBold(outputAttributes,true);
							}
						if(b2==true)
							{
							StyleConstants.setItalic(outputAttributes,true);
							}
						StyleConstants.setFontSize(outputAttributes,s);
						outputDocument.insertString(outputDocument.getLength(),stri,outputAttributes);
//						repaint();
						StyleConstants.setBold(outputAttributes,false);
						StyleConstants.setItalic(outputAttributes,false);
						StyleConstants.setFontSize(outputAttributes,scriptsize);
						StyleConstants.setForeground(outputAttributes,Color.black);
						}
					catch(BadLocationException ble)
						{
						System.err.print(ble.getMessage());
						}
					}
				});
			}
		catch(InterruptedException ie)
			{
			System.err.print(ie.getMessage());
			}
		catch(java.lang.reflect.InvocationTargetException ite)
			{
			System.err.print(ite.getMessage());
			}
		}
	}

private void JavaOutput(String str)
        {
	if(consoleInput)
		{
		System.out.print(str);
		}
	else
		{
		try
			{
			final String stri = str;
			SwingUtilities.invokeAndWait(new Runnable()
				{
				public void run()
					{
					try
						{
						outputDocument.insertString(outputDocument.getLength(),stri,outputAttributes);
//						repaint();
						}
					catch(BadLocationException ble)
						{
						System.err.print(ble.getMessage());
						}
					}
				});
			}
		catch(InterruptedException ie)
			{
			System.err.print(ie.getMessage());
			}
		catch(java.lang.reflect.InvocationTargetException ite)
			{
			System.err.print(ite.getMessage());
			}
		}
	}


public void Out(String str, boolean thick, boolean italic, short size, int c1, int c2, int c3)
        {
	if(consoleInput)
		{
		System.out.print(str);
		}
	else
		{
		final String stri = str;
		final boolean b1 = thick;
		final boolean b2 = italic;
		final int s = (int)(size+fontFactor);
		final int i1 = c1;
		final int i2 = c2;
		final int i3 = c3;
		try
			{
			if(i1!=0 || i2!=0 || i3!=0)
				{
				StyleConstants.setForeground(outputAttributes, new Color(i1, i2, i3));
				}
			if(b1==true)
				{
				StyleConstants.setBold(outputAttributes,true);
				}
			if(b2==true)
				{
				StyleConstants.setItalic(outputAttributes,true);
				}
			StyleConstants.setFontSize(outputAttributes,s);
			outputDocument.insertString(outputDocument.getLength(),stri,outputAttributes);
			StyleConstants.setBold(outputAttributes,false);
			StyleConstants.setItalic(outputAttributes,false);
			StyleConstants.setFontSize(outputAttributes,scriptsize);
			StyleConstants.setForeground(outputAttributes,Color.black);
			}
		catch(BadLocationException ble)
			{
			System.err.print(ble.getMessage());
			}
		}
	}


public void Out(String str)
        {
	if(consoleInput)
		{
		System.out.print(str);
		}
	else
		{
		final String stri = str;
		try
			{
			outputDocument.insertString(outputDocument.getLength(),stri,outputAttributes);
			}
		catch(BadLocationException ble)
			{
			System.err.print(ble.getMessage());
			}
		}
	}

private void JavaReview(String str)
	{
	if(!consoleInput)
		{
		if(reviewVector.size()>100)
			{
			reviewVector.removeElementAt(0);
			}
		reviewVector.add(str);
		reviewList.setListData(reviewVector);
		}
	}

private void ClearOutput()
	{
	fileCommand(2);
	}


private void addtoVector(Vector v, String str)
	{
	v.add(str);
	}

private void JavaSaveOutput()
	{
	saveOutput();
	}


private void setDelim(String str)
        {
        if (str.equals("return"))
          {
          delimiter = '\n';
          }
        else if (str.equals("semicolon"))
          {
          delimiter = ';';
          }
        }


private native void setObjectList(Vector v, String type);
private native void setObjectTypeList(Vector v);

// fuer describe dataset
private native String getValue(int i, int j);
protected native void setVarnames(Vector v);
protected native int getRows();

// fuer graphobj
protected native double getDoubleValue(int i, int j);
protected native int getDRows();
protected native int getDCols();
protected native double getMax(int col);
protected native double getMin(int col);
protected native String getVarname(int nr);

// buttons
private native void setStop(boolean stop);
protected native void setPause(boolean pause);
native void setProcessrunning(boolean proc);
private native void setSuppressoutput(boolean supp);

// random number generation

private double juniform()
	{
	return UniformGen.nextDouble(rStream, 0, 1);
	}

private double jnormal()
	{
	return NormalGen.nextDouble(rStream, 0, 1);
	}

private double jexponential(double rate)
	{
	return ExponentialGen.nextDouble(rStream, rate);
	}

private double jbernoulli(double pi)
	{
	if(pi==1)
		return 1;
	else
		return BinomialGen.nextInt(rStream, 1, pi);
	}

private double jbinomial(double size, double pi)
	{
	if(pi==1)
		return (int)size;
	else
		return BinomialGen.nextInt(rStream, (int)size, pi);
	}

private double jgamma(double mu, double nu)
	{
	return GammaGen.nextDouble(rStream, nu, nu/mu);
	}

private double jpoisson(double lambda)
	{
	return PoissonGen.nextInt(rStream, lambda);
	}

private double jweibull(double alpha, double lambda)
	{
	return WeibullGen.nextDouble(rStream, alpha, lambda, 0.0);
	}


//Funktion, die bei der uebergabe eines Befehls aufgerufen wird

private void doparse(String inp)
	{
	breakButton.setEnabled(true);
	setStop(false);
	do
		{
		int i = inp.indexOf(delimiter);
		inp = inp.substring(0,Math.max(i,0))+inp.substring(Math.min(i+1,inp.length()),inp.length());
		}
		while(inp.indexOf(delimiter)>-1);
	commandArea.setText("");
	if(!inp.equals("quit") && !inp.equals("exit"))
		{
		if(inp.length()>0 && inp.charAt(0)!='%')
			{
			Out("> "+inp+"\n");
			if(reviewVector.size()>100)
				{
				reviewVector.removeElementAt(0);
				}
			reviewVector.add(inp);
			reviewList.setListData(reviewVector);
			}
		}
	setProcessrunning(true);
	processRunning=true;
	commandArea.setEnabled(false);
	commandArea.getCaret().setVisible(false);
	comm.setEnabled(false);
	objectsList1.setEnabled(false);
	objectsList2.setEnabled(false);
	obj.setEnabled(false);
	open.setEnabled(false);
	parseThread = new ParseThread(this);
	parseThread.setPriority(threadPriority);
	parseThread.setCommand(inp);
	outputPane.setCaretPosition(outputDocument.getLength());
	if(pause)
		{
		Out("\nPROGRAM PAUSED\nClick CONTINUE to proceed\n\n");
		parseThread.start();
		parseThread.setPriority(Thread.MIN_PRIORITY);
//		parseThread.suspend();
		}
	else
		{
		parseThread.start();
		}
	}

//Hilfsfunktionen die ParseThread verwendet

void finished()
	{
	processRunning=false;
	commandArea.setEnabled(true);
	comm.setEnabled(true);
	objectsList1.setEnabled(true);
	objectsList2.setEnabled(true);
	obj.setEnabled(true);
	open.setEnabled(true);
	reviewScrollPane.getVerticalScrollBar().setValue(reviewScrollPane.getVerticalScrollBar().getMaximum()+reviewScrollPane.getVerticalScrollBar().getVisibleAmount());
	reviewScrollPane.getHorizontalScrollBar().setValue(0);
	commandArea.getCaret().setVisible(true);
	breakButton.setEnabled(false);
	commandArea.setEnabled(true);
	commandArea.getCaret().setVisible(true);
	commandArea.requestFocus();
	int index1 = objectsList1.getSelectedIndex();
	if(index1>-1)
		{
		objectVector.clear();
		String selected = (String)objectsVector.elementAt(index1);
		setObjectList(objectVector,selected);
		objectsList2.setListData(objectVector);
		}
	}


native boolean parse(String str);
private native void parsecommand(String str);
private native void parsecommand2(JTextArea oA, String str);


//Klasse main, die ein BayesX-Fenster erzeugt und anzeigt
public static void main(String[] args)
	{
	if(args.length==1 && args[0].equals("-g"))
		{
		BayesX b = new BayesX();
		b.consoleInput = true;
		System.out.println("\nBayesX Version 2.1");
		System.out.println("Software for Bayesian Inference in Structured Additive Regression Models");
		BufferedReader console = new BufferedReader(new InputStreamReader(System.in));
		String str = null;
		while(true)
			{
			try
				{
				System.out.print("\n> ");
				str = console.readLine();
				if(str.equals("exit")||str.equals("quit"))
					{
					System.exit(0);
					}
				else
					{
					b.parsecommand(str);
					}
				}
			catch(IOException ioe)
				{
				System.err.println(ioe.getMessage());
				}
			}
		}
	else
		{
	 	BayesX b = new BayesX();
		b.show();
		}
	}


}




