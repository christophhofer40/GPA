/*
this plugin uses the geometric phase analysis to map the strain of HR-(S)TEM images
also see:
Hytch, Ultramicroscopy 1998
*/

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import ij.plugin.filter.*;
import static ij.measure.Measurements.MEAN;

public class GPA implements PlugIn
{
	//metadata of input image
	ImagePlus imp;
	ImageProcessor impr;
	FHT fht;	
	int maxN;
	int originalWidth;
	int originalHeight;
	boolean padded;

	



	public void run(String arg) 
	{
		imp = IJ.getImage();
		impr = imp.getProcessor();
		FHT fht=newFHT(impr);
		fht.transform();
		ImageStack complexSt = fht.getComplexTransform();
		ImagePlus complexP = new ImagePlus("Complex Stack",complexSt);
		ImagePlus PS=new ImagePlus("PS of Image",fht.getPowerSpectrum());
		PS.show();
		if(!doDialog()) return;
		Roi roiPS=PS.getRoi();
		Roi roiIm=imp.getRoi();
		while(roiIm==null ||roiPS==null || 
			roiPS.getType()!=Roi.OVAL|| roiIm.getType()!=Roi.RECTANGLE)
		{
			roiPS=PS.getRoi();
			roiIm=imp.getRoi();
			if(!doDialog()) return;
		}
		
		FHT refFHT = new FHT(pad(imp.duplicate().getProcessor(),maxN));
		refFHT.transform();
		ImagePlus refPS=new ImagePlus("PS of reference",refFHT.getPowerSpectrum());
		refPS.setRoi(roiPS);
		ImageStatistics stats = ImageStatistics.getStatistics(refPS.getProcessor());
		System.out.println(stats.max);
		refPS.show();
		//IJ.doCommand("Make Inverse");
		/*


		float[] re=(float[])imp.getStack().getPixels(1);
		float[] im=(float[])imp.getStack().getPixels(2);
		float[] phase=new float[re.length];
		int imsize=re.length; 
		for (int i=0;i<imsize;++i)
		{	
			float ph=(float)Math.atan2(im[i],re[i]);
			phase[i]=ph;
		}
		ImageStack imst=new ImageStack(imp.getWidth(),imp.getHeight(),1);
		imst.setPixels(phase,1);
		ImagePlus imphase=new ImagePlus("Phase",imst);
		imphase.show();

		*/
	}

	boolean doDialog()
	{
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(getClass().getSimpleName());
		gd.addMessage("Please specify aperture in the Power Spectrum and a reference in the real space Image");
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		return true;
	}


	ImageProcessor pad(ImageProcessor ip) 
    	{
        		originalWidth = ip.getWidth();
      		originalHeight = ip.getHeight();
		maxN = Math.max(originalWidth, originalHeight);
        		int i = 2;
        		while(i<maxN) i *= 2;
        		if (i==maxN && originalWidth==originalHeight) 
		{
            			padded = false;
           		 	return ip;
        		}
        		maxN = i;
        		//showStatus("Padding to "+ maxN + "x" + maxN);
        		ImageStatistics stats = ImageStatistics.getStatistics(ip, MEAN, null);
        		ImageProcessor ip2 = ip.createProcessor(maxN, maxN);
        		ip2.setValue(stats.mean);
        		ip2.fill();
        		ip2.insert(ip, 0, 0);
        		padded = true;
        		Undo.reset();
        		//new ImagePlus("padded", ip2.duplicate()).show();
        		return ip2;
    	}
	
	//pad with a certain size
	ImageProcessor pad(ImageProcessor ip, int maxN) 
    	{
        		//showStatus("Padding to "+ maxN + "x" + maxN);
        		ImageStatistics stats = ImageStatistics.getStatistics(ip, MEAN, null);
        		ImageProcessor ip2 = ip.createProcessor(maxN, maxN);
        		ip2.setValue(stats.mean);
        		ip2.fill();
        		ip2.insert(ip, 0, 0);
        		padded = true;
        		Undo.reset();
        		//new ImagePlus("padded", ip2.duplicate()).show();
        		return ip2;
    	}



	FHT newFHT(ImageProcessor ip) 
   	{
       		FHT fht;
       		if (ip instanceof ColorProcessor) 
        		{
            			//showStatus("Extracting brightness");
            			ImageProcessor ip2 = ((ColorProcessor)ip).getBrightness();
           			fht = new FHT(pad(ip2));
            			fht.rgb = (ColorProcessor)ip.duplicate(); // save so we can later update the brightness
        		}	 
       		else fht = new FHT(pad(ip));
       		if (padded) 
		{
           			fht.originalWidth = originalWidth;
           			fht.originalHeight = originalHeight;
        		}	
       		int bitDepth = imp.getBitDepth();
        		fht.originalBitDepth = bitDepth;
        		if (bitDepth!=24)
        			fht.originalColorModel = ip.getColorModel();
       	 	return fht;
    	}
    
   	float[] crop(float[] ar1, int w1, int w2, int h2)
    	{
      		float[] ar2= new float[w2*h2];
        		//int h1=ar1.length/w1;
       		int c=0;
        		for(int i=0;i<ar1.length;++i)
       		{
           			int x=(i%w1);
           			int y=(int)(i/w1);
           			if(x>=w2 || y>=h2) continue;
            			ar2[c]=ar1[i];
            			c++;
       	 	}
        		return ar2;
	}


}
