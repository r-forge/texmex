
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body style="background-color:FFFF99; color:0000FF">

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<h1 align="center"><tt>texmex</tt></h1>

<h1><u><tt>t</tt></u>hreshold <u><tt>ex</tt></u>ceedences and <u><tt>m</tt></u>ultivariate <u><tt>ex</tt></u>tremes</h1>

<p> This project will produce an R package, <tt>texmex</tt>,
for modelling extreme values
using generalized Pareto distributions. Inlcuded will be maximum likelihood
and penalized likelihood (maximum a posteriori) estimation for threshold
exceedences, including a formula interface. Also included will be a fully
Bayesian estimation algorithm using MCMC. Model objects will have methods
for plotting and summarizing the models. </p>

<!-- <table>
  <tbody>
    <tr>
      <td> -->
          <img src="bwinter.png" alt="Posterior distributions and diagnostic plots"/>
<!--      </td>
      <td> -->

<!--      </td>
    </tr>
  </tbody>
</table> -->

<p> Also included will be an implementation of the Heffernan-Tawn (RSS B, 2004)
approach to conditional multivariate extreme value modelling, again with plot
and other methods.</p>

<p>We intend to include a comprehensive suite of test scripts.</p>

          <img src="pWinter.png" alt="Predictions from a conditional multivariate extreme value model"/>

<p>The code currently included has been written by Janet E. Heffernan and
Harry Southworth, with additional material by Ioannis Papastathopoulos.
The <tt>chi</tt> functions are wrapped versions of functions appearing
in the <tt>evd</tt> package, maintained by Alec Stephenson.</p>

<p><strong>There are several known major bugs and problems.</strong>
This project is currently closed to new developers. Once we get a version
out of beta, please feel free to ask to join.</p>




<h3>References</h3>
<p>J. E. Heffernan and J. A. Tawn, A conditional approach to multivariate
extreme values, Journal of the Royal Statistical Society (B), 66, 
497 - 546, 2004</p>

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>
