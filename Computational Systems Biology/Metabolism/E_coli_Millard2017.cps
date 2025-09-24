<?xml version="1.0" encoding="UTF-8"?>
<!-- generated with COPASI 4.15 (Build 95) (http://www.copasi.org) at 2016-05-31 21:44:46 UTC -->
<?oxygen RNGSchema="http://www.copasi.org/static/schema/CopasiML.rng" type="xml"?>
<COPASI xmlns="http://www.copasi.org/static/schema" versionMajor="4" versionMinor="15" versionDevel="95" copasiSourcesModified="0">
  <ListOfFunctions>
    <Function key="Function_13" name="Mass action (irreversible)" type="MassAction" reversible="false">
      <MiriamAnnotation>
<rdf:RDF xmlns:CopasiMT="http://www.copasi.org/RDF/MiriamTerms#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
   <rdf:Description rdf:about="#Function_13">
   <CopasiMT:is rdf:resource="urn:miriam:obo.sbo:SBO:0000041" />
   </rdf:Description>
   </rdf:RDF>
      </MiriamAnnotation>
      <Comment>
        <body xmlns="http://www.w3.org/1999/xhtml">
<b>Mass action rate law for first order irreversible reactions</b>
<p>
Reaction scheme where the products are created from the reactants and the change of a product quantity is proportional to the product of reactant activities. The reaction scheme does not include any reverse process that creates the reactants from the products. The change of a product quantity is proportional to the quantity of one reactant.
</p>
</body>
      </Comment>
      <Expression>
        k1*PRODUCT&lt;substrate_i&gt;
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_81" name="k1" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_79" name="substrate" order="1" role="substrate"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_14" name="Mass action (reversible)" type="MassAction" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:CopasiMT="http://www.copasi.org/RDF/MiriamTerms#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
   <rdf:Description rdf:about="#Function_14">
   <CopasiMT:is rdf:resource="urn:miriam:obo.sbo:SBO:0000042" />
   </rdf:Description>
   </rdf:RDF>
      </MiriamAnnotation>
      <Comment>
        <body xmlns="http://www.w3.org/1999/xhtml">
<b>Mass action rate law for reversible reactions</b>
<p>
Reaction scheme where the products are created from the reactants and the change of a product quantity is proportional to the product of reactant activities. The reaction scheme does include a reverse process that creates the reactants from the products.
</p>
</body>
      </Comment>
      <Expression>
        k1*PRODUCT&lt;substrate_i&gt;-k2*PRODUCT&lt;product_j&gt;
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_62" name="k1" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_39" name="substrate" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_67" name="k2" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_71" name="product" order="3" role="product"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_40" name="Function for PGI" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_40">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-05T14:27:25Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(G6P-F6P/Keq)/(KmG6P)/(1+F6P/KmF6P+G6P/KmG6P+PEP/KmPEP+PGN/KmPGN)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_258" name="F6P" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_254" name="G6P" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_266" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_270" name="KmF6P" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_268" name="KmG6P" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_274" name="KmPEP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_272" name="PEP" order="6" role="modifier"/>
        <ParameterDescription key="FunctionParameter_264" name="Vmax" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_276" name="PGN" order="8" role="modifier"/>
        <ParameterDescription key="FunctionParameter_278" name="KmPGN" order="9" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_41" name="Function for PFK" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_41">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2014-12-29T22:03:23Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*n*(MgATP*F6P-MgADP*FDP/Keq)/(KirF6P*KmrATPMg)/(1+KmrFDP/KirFDP*(MgADP/KmrADP)+KmrF6P/KirF6P*(MgATP/KmrATPMg)+KmrFDP/KirFDP*(MgADP/KmrADP)*(F6P/KirF6P)+MgATP/KmrATPMg*(F6P/KirF6P)+MgADP/KirADP*(MgATP/KmrATPMg)*(F6P/KirF6P)+(1+(ATP-MgATP)/KirATP)*(F6P/KirF6P)+FDP/KirFDP+MgADP/KmrADP*(FDP/KirFDP)+KmrF6P/KirF6P*(MgATP/KmrATPMg)*(FDP/KirFDP)+Wr*(KmrF6P/KirF6P)*(MgADP/KirADP)*(MgATP/KmrATPMg)*(FDP/KmrFDP))/(1+L0*((1+KmtFDP/KitFDP*(MgADP/KmtADP)+KmtF6P/KitF6P*(MgATP/KmtATPMg)+KmtFDP/KitFDP*(MgADP/KmtADP)*(F6P/KitF6P)+MgATP/KmtATPMg*(F6P/KitF6P)+MgADP/KitADP*(MgATP/KmtATPMg)*(F6P/KitF6P)+(1+(ATP-MgATP)/KitATP)*(F6P/KitF6P)+FDP/KitFDP+MgADP/KmtADP*(FDP/KitFDP)+KmtF6P/KitF6P*(MgATP/KmtATPMg)*(FDP/KitFDP)+Wt*(KmtF6P/KitF6P)*(MgADP/KitADP)*(MgATP/KmtATPMg)*(FDP/KmtFDP))*(1+MgADP/KeftADP+PEP/KeftPEP+MgADP/KeftADP*(PEP/KeftPEP))/((1+KmrFDP/KirFDP*(MgADP/KmrADP)+KmrF6P*MgATP/(KirF6P*KmrATPMg)+KmrFDP/KirFDP*(MgADP/KmrADP)*(F6P/KirF6P)+MgATP/KmrATPMg*(F6P/KirF6P)+MgADP/KirADP*(MgATP/KmrATPMg)*(F6P/KirF6P)+(1+(ATP-MgATP)/KirATP)*(F6P/KirF6P)+FDP/KirFDP+MgADP/KmrADP*(FDP/KirFDP)+KmrF6P/KirF6P*(MgATP/KmrATPMg)*(FDP/KirFDP)+Wr*(KmrF6P/KirF6P)*(MgADP/KirADP)*(MgATP/KmrATPMg)*(FDP/KmrFDP))*(1+MgADP/KefrADP+PEP/KefrPEP+MgADP/KefrADP*(PEP/KefrPEP))))^n)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_298" name="ATP" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_267" name="F6P" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_280" name="FDP" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_330" name="KefrADP" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_332" name="KefrPEP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_324" name="KeftADP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_328" name="KeftPEP" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_282" name="Keq" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_296" name="KirADP" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_300" name="KirATP" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_284" name="KirF6P" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_290" name="KirFDP" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_318" name="KitADP" order="12" role="constant"/>
        <ParameterDescription key="FunctionParameter_320" name="KitATP" order="13" role="constant"/>
        <ParameterDescription key="FunctionParameter_314" name="KitF6P" order="14" role="constant"/>
        <ParameterDescription key="FunctionParameter_308" name="KitFDP" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_292" name="KmrADP" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_286" name="KmrATPMg" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_294" name="KmrF6P" order="18" role="constant"/>
        <ParameterDescription key="FunctionParameter_288" name="KmrFDP" order="19" role="constant"/>
        <ParameterDescription key="FunctionParameter_310" name="KmtADP" order="20" role="constant"/>
        <ParameterDescription key="FunctionParameter_316" name="KmtATPMg" order="21" role="constant"/>
        <ParameterDescription key="FunctionParameter_312" name="KmtF6P" order="22" role="constant"/>
        <ParameterDescription key="FunctionParameter_306" name="KmtFDP" order="23" role="constant"/>
        <ParameterDescription key="FunctionParameter_304" name="L0" order="24" role="constant"/>
        <ParameterDescription key="FunctionParameter_246" name="MgADP" order="25" role="modifier"/>
        <ParameterDescription key="FunctionParameter_271" name="MgATP" order="26" role="modifier"/>
        <ParameterDescription key="FunctionParameter_326" name="PEP" order="27" role="modifier"/>
        <ParameterDescription key="FunctionParameter_279" name="Vmax" order="28" role="constant"/>
        <ParameterDescription key="FunctionParameter_302" name="Wr" order="29" role="constant"/>
        <ParameterDescription key="FunctionParameter_322" name="Wt" order="30" role="constant"/>
        <ParameterDescription key="FunctionParameter_275" name="n" order="31" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_42" name="Function for FBA" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_42">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2014-12-29T22:00:27Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(FDP-DAP*GAP/Keq)/KmFDP/(1+FDP/KmFDP+DAP/KmDAP+DAP/KmDAP*(GAP/KmGAP)+PEP/KmPEP)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_325" name="DAP" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_329" name="FDP" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_321" name="GAP" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_317" name="Keq" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_309" name="KmDAP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_313" name="KmFDP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_305" name="KmGAP" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_297" name="KmPEP" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_301" name="PEP" order="8" role="modifier"/>
        <ParameterDescription key="FunctionParameter_333" name="Vmax" order="9" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_43" name="Function for TPI" type="UserDefined" reversible="true">
      <Expression>
        Vmax*(DAP-GAP/Keq)/KmDAP/(1+DAP/KmDAP+GAP/KmGAP)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_303" name="DAP" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_311" name="GAP" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_319" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_327" name="KmDAP" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_293" name="KmGAP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_295" name="Vmax" order="5" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_44" name="Function for GDH" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_44">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:43:08Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(P*GAP*NAD-BPG*NADH/Keq)/(KmP*KmGAP*KmNAD)/((1+P/KmP)*(1+GAP/KmGAP)*(1+NAD/KmNAD)+(1+BPG/KmBPG)*(1+NADH/KmNADH)-1)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_285" name="BPG" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_307" name="GAP" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_265" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_338" name="KmBPG" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_334" name="KmGAP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_336" name="KmNAD" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_340" name="KmNADH" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_273" name="KmP" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_289" name="NAD" order="8" role="substrate"/>
        <ParameterDescription key="FunctionParameter_281" name="NADH" order="9" role="product"/>
        <ParameterDescription key="FunctionParameter_323" name="P" order="10" role="substrate"/>
        <ParameterDescription key="FunctionParameter_291" name="Vmax" order="11" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_45" name="Function for PGK" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_45">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T08:49:30Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(MgADP*BPG-MgATP*PGA3/Keq)/(KmADPMg*KmBPG)/(1+MgADP/KmADPMg+BPG/KmBPG+MgADP/KmADPMg*BPG/KmBPG+MgATP/KmATPMg+PGA3/KmPGA3+MgATP/KmATPMg*PGA3/KmPGA3)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_277" name="BPG" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_315" name="Keq" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_342" name="KmADPMg" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_346" name="KmATPMg" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_344" name="KmBPG" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_348" name="KmPGA3" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_337" name="MgADP" order="6" role="modifier"/>
        <ParameterDescription key="FunctionParameter_262" name="MgATP" order="7" role="modifier"/>
        <ParameterDescription key="FunctionParameter_287" name="PGA3" order="8" role="product"/>
        <ParameterDescription key="FunctionParameter_341" name="Vmax" order="9" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_46" name="Function for GPM" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_46">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-13T10:20:43Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(PGA3-PGA2/Keq)/KmPGA3/(1+PGA3/KmPGA3+PGA2/KmPGA2)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_283" name="Keq" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_350" name="KmPGA2" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_335" name="KmPGA3" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_331" name="PGA2" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_345" name="PGA3" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_349" name="Vmax" order="5" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_47" name="Function for ENO" type="UserDefined" reversible="true">
      <Expression>
        Vmax*(PGA2-PEP/Keq)/KmPGA2/(1+PGA2/KmPGA2+PEP/KmPEP)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_352" name="Keq" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_356" name="KmPEP" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_354" name="KmPGA2" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_343" name="PEP" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_269" name="PGA2" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_351" name="Vmax" order="5" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_48" name="Function for PYK" type="UserDefined" reversible="false">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_48">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T08:51:41Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*n*PEP*MgADP/(KirPEP*KmrADPMg)/(1+KmrPEP/KirPEP*(MgADP/KmrADPMg)+MgATP/KirATP+MgADP/KmrADPMg*(PEP/KirPEP)+KmrADPMg/KmrADPMg*(1+(ADP-MgADP)/KirADP)*(PEP/KirPEP)+PYR/KirPYR+MgATP/KirPyrATP*(PYR/KirPYR))/(1+L0*((1+KmtPEP/KitPEP*(MgADP/KmtADPMg)+MgATP/KitATP+MgADP*PEP/(KitPEP*KmtADPMg)+(1+(ADP-MgADP)/KitADP)*(PEP/KitPEP)+PYR/KitPYR+MgATP/KitPyrATP*(PYR/KitPYR))*(1+SUCCOA/KeftSUCCOA+MgATP*SUCCOA/(KeftATP*KeftSUCCOA))/((1+KmrPEP/KirPEP*(MgADP/KmrADPMg)+MgATP/KirATP+MgADP/KmrADPMg*(PEP/KirPEP)+(1+(ADP-MgADP)/KirADP)*(PEP/KirPEP)+PYR/KirPYR+MgATP/KirPyrATP*(PYR/KirPYR))*(1+FDP/KefrFDP+G6P/KefrG6P+GL6P/KefrGL6P+R5P/KefrR5P+RU5P/KefrRU5P+S7P/KefrS7P+X5P/KefrX5P)))^n)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_370" name="ADP" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_402" name="FDP" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_406" name="G6P" order="2" role="modifier"/>
        <ParameterDescription key="FunctionParameter_410" name="GL6P" order="3" role="modifier"/>
        <ParameterDescription key="FunctionParameter_404" name="KefrFDP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_408" name="KefrG6P" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_412" name="KefrGL6P" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_416" name="KefrR5P" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_420" name="KefrRU5P" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_424" name="KefrS7P" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_428" name="KefrX5P" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_400" name="KeftATP" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_398" name="KeftSUCCOA" order="12" role="constant"/>
        <ParameterDescription key="FunctionParameter_372" name="KirADP" order="13" role="constant"/>
        <ParameterDescription key="FunctionParameter_368" name="KirATP" order="14" role="constant"/>
        <ParameterDescription key="FunctionParameter_360" name="KirPEP" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_376" name="KirPYR" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_378" name="KirPyrATP" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_390" name="KitADP" order="18" role="constant"/>
        <ParameterDescription key="FunctionParameter_388" name="KitATP" order="19" role="constant"/>
        <ParameterDescription key="FunctionParameter_384" name="KitPEP" order="20" role="constant"/>
        <ParameterDescription key="FunctionParameter_392" name="KitPYR" order="21" role="constant"/>
        <ParameterDescription key="FunctionParameter_394" name="KitPyrATP" order="22" role="constant"/>
        <ParameterDescription key="FunctionParameter_362" name="KmrADPMg" order="23" role="constant"/>
        <ParameterDescription key="FunctionParameter_364" name="KmrPEP" order="24" role="constant"/>
        <ParameterDescription key="FunctionParameter_386" name="KmtADPMg" order="25" role="constant"/>
        <ParameterDescription key="FunctionParameter_382" name="KmtPEP" order="26" role="constant"/>
        <ParameterDescription key="FunctionParameter_380" name="L0" order="27" role="constant"/>
        <ParameterDescription key="FunctionParameter_358" name="MgADP" order="28" role="modifier"/>
        <ParameterDescription key="FunctionParameter_366" name="MgATP" order="29" role="modifier"/>
        <ParameterDescription key="FunctionParameter_299" name="PEP" order="30" role="substrate"/>
        <ParameterDescription key="FunctionParameter_374" name="PYR" order="31" role="product"/>
        <ParameterDescription key="FunctionParameter_414" name="R5P" order="32" role="modifier"/>
        <ParameterDescription key="FunctionParameter_418" name="RU5P" order="33" role="modifier"/>
        <ParameterDescription key="FunctionParameter_422" name="S7P" order="34" role="modifier"/>
        <ParameterDescription key="FunctionParameter_396" name="SUCCOA" order="35" role="modifier"/>
        <ParameterDescription key="FunctionParameter_357" name="Vmax" order="36" role="constant"/>
        <ParameterDescription key="FunctionParameter_426" name="X5P" order="37" role="modifier"/>
        <ParameterDescription key="FunctionParameter_353" name="n" order="38" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_49" name="Function for ZWF" type="UserDefined" reversible="true">
      <Expression>
        Vmax*(G6P*NADP-GL6P*NADPH/Keq)/(KdG6P*KmNADP)/(1+G6P/KdG6P+KmG6P/KdG6P*(NADP/KmNADP)+G6P/KdG6P*(NADP/KmNADP)+KmGL6P/KdGL6P*(NADPH/KmNADPH)+GL6P/KdGL6P*(NADPH/KmNADPH))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_425" name="G6P" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_417" name="GL6P" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_405" name="KdG6P" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_389" name="KdGL6P" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_409" name="Keq" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_397" name="KmG6P" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_393" name="KmGL6P" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_401" name="KmNADP" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_385" name="KmNADPH" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_421" name="NADP" order="9" role="substrate"/>
        <ParameterDescription key="FunctionParameter_413" name="NADPH" order="10" role="product"/>
        <ParameterDescription key="FunctionParameter_429" name="Vmax" order="11" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_50" name="Function for PGL" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_50">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2014-12-29T22:03:27Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(GL6P-PGN/Keq)/KmGL6P/(1+GL6P/KmGL6P+PGN/KmPGN+G6P/KiG6P)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_381" name="G6P" order="0" role="modifier"/>
        <ParameterDescription key="FunctionParameter_391" name="GL6P" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_407" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_377" name="KiG6P" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_415" name="KmGL6P" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_423" name="KmPGN" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_399" name="PGN" order="6" role="product"/>
        <ParameterDescription key="FunctionParameter_383" name="Vmax" order="7" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_51" name="Function for GND" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_51">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:50:19Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(NADP*PGN-NADPH*RU5P*HCO3/Keq)/(KdNADP*KmPGN)/(1+NADP/KdNADP+FDP/KefFbP+NADP/KdNADP*(FDP/KefNADPFbP)+KmNADP/KdNADP*(PGN/KmPGN)+NADP/KdNADP*(PGN/KmPGN)+ATP/KefATP+ATP/KefNADPATP*(KmNADP/KdNADP)*(PGN/KmPGN)+HCO3/KdHCO3+NADPH/KdNADPH+RU5P/KdRu5P+HCO3/KdHCO3*(NADPH/KdHCO3NADPH)+HCO3/KdHCO3*(KmNADPH/KdHCO3NADPH)*(RU5P/KmRU5P)+HCO3/KdHCO3*(NADPH/KdHCO3NADPH)*(RU5P/KmRU5P)+KmHCO3/KdHCO3*(NADPH/KdHCO3NADPH*(RU5P/KmRU5P)))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_437" name="ATP" order="0" role="modifier"/>
        <ParameterDescription key="FunctionParameter_355" name="FDP" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_369" name="HCO3" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_443" name="KdHCO3" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_449" name="KdHCO3NADPH" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_361" name="KdNADP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_445" name="KdNADPH" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_447" name="KdRu5P" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_439" name="KefATP" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_431" name="KefFbP" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_441" name="KefNADPATP" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_433" name="KefNADPFbP" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_365" name="Keq" order="12" role="constant"/>
        <ParameterDescription key="FunctionParameter_455" name="KmHCO3" order="13" role="constant"/>
        <ParameterDescription key="FunctionParameter_435" name="KmNADP" order="14" role="constant"/>
        <ParameterDescription key="FunctionParameter_451" name="KmNADPH" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_339" name="KmPGN" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_453" name="KmRU5P" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_427" name="NADP" order="18" role="substrate"/>
        <ParameterDescription key="FunctionParameter_395" name="NADPH" order="19" role="product"/>
        <ParameterDescription key="FunctionParameter_411" name="PGN" order="20" role="substrate"/>
        <ParameterDescription key="FunctionParameter_373" name="RU5P" order="21" role="product"/>
        <ParameterDescription key="FunctionParameter_375" name="Vmax" order="22" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_52" name="Function for RPE" type="UserDefined" reversible="true">
      <Expression>
        Vmax*(RU5P-X5P/Keq)/KmRU5P/(1+RU5P/KmRU5P+X5P/KmX5P)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_444" name="Keq" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_440" name="KmRU5P" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_436" name="KmX5P" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_452" name="RU5P" order="3" role="substrate"/>
        <ParameterDescription key="FunctionParameter_456" name="Vmax" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_448" name="X5P" order="5" role="product"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_53" name="Function for RPI" type="UserDefined" reversible="true">
      <Expression>
        Vmax*(RU5P-R5P/Keq)/KmRU5P/(1+RU5P/KmRU5P+R5P/KmR5P+E4P/KmE4P)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_371" name="E4P" order="0" role="modifier"/>
        <ParameterDescription key="FunctionParameter_432" name="Keq" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_403" name="KmE4P" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_363" name="KmR5P" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_347" name="KmRU5P" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_450" name="R5P" order="5" role="product"/>
        <ParameterDescription key="FunctionParameter_442" name="RU5P" order="6" role="substrate"/>
        <ParameterDescription key="FunctionParameter_434" name="Vmax" order="7" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_54" name="Function for X5P_GAP_TKT" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_54">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-31T12:39:06Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        kcat*(tkt*X5P-GAP*tktC2/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_446" name="GAP" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_458" name="Keq" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_430" name="X5P" order="2" role="substrate"/>
        <ParameterDescription key="FunctionParameter_419" name="kcat" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_367" name="tkt" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_379" name="tktC2" order="5" role="product"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_55" name="Function for F6P_E4P_TKT" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_55">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2014-12-29T22:00:24Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        kcat*(E4P*tktC2-F6P*tkt/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_438" name="E4P" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_460" name="F6P" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_464" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_459" name="kcat" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_462" name="tkt" order="4" role="product"/>
        <ParameterDescription key="FunctionParameter_359" name="tktC2" order="5" role="substrate"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_56" name="Function for S7P_R5P_TKT" type="UserDefined" reversible="true">
      <Expression>
        kcat*(R5P*tktC2-S7P*tkt/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_470" name="Keq" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_461" name="R5P" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_466" name="S7P" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_465" name="kcat" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_468" name="tkt" order="4" role="product"/>
        <ParameterDescription key="FunctionParameter_454" name="tktC2" order="5" role="substrate"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_57" name="Function for F6P_GAP_TAL" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_57">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2014-12-29T22:00:25Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        kcat*(GAP*talC3-F6P*tal/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_472" name="F6P" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_467" name="GAP" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_476" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_471" name="kcat" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_474" name="tal" order="4" role="product"/>
        <ParameterDescription key="FunctionParameter_387" name="talC3" order="5" role="substrate"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_58" name="Function for S7P_E4P_TAL" type="UserDefined" reversible="true">
      <Expression>
        kcat*(S7P*tal-E4P*talC3/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_478" name="E4P" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_482" name="Keq" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_473" name="S7P" order="2" role="substrate"/>
        <ParameterDescription key="FunctionParameter_477" name="kcat" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_457" name="tal" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_480" name="talC3" order="5" role="product"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_59" name="Function for FBP" type="UserDefined" reversible="false">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_59">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:42:50Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*n*MgFDP/KirFDPMg/(1+KmrFDP/KirFDP*(MG/KmrMg)+P/KirP+P/KirP*(MG/KirPMg)+F6P/KirF6P+F6P/KirF6P*(MG/KirF6PMg)+P/KirP*(F6P/KirPF6P)+P/KirP*(F6P/KirPF6P)*(MG/KirPF6PMg)+(FDP-MgFDP)/KirFDP+KdFDPMg/KmrMg*(MgFDP/KirFDP)+AMP/KirAMP+MgFDP/KirFDPMg+MgFDP/KirFDPMg*(MG/KirFDPMgMg)+AMP/KirAMP*((FDP-MgFDP)/KirAMPFDP))/(1+L0*((1+KmtFDP/KitFDP*(MG/KmtMg)+P/KitP+P/KitP*(MG/KitPMg)+F6P/KitF6P+F6P/KitF6P*(MG/KitF6PMg)+P/KitP*(F6P/KitPF6P)+P/KitP*(F6P/KitPF6P)*(MG/KitPF6PMg)+(FDP-MgFDP)/KitFDP+KdFDPMg/KmtMg*(MgFDP/KitFDP)+AMP/KitAMP+MgFDP/KitFDPMg+MgFDP/KitFDPMg*(MG/KitFDPMgMg)+AMP/KitAMP*((FDP-MgFDP)/KitAMPFDP))/(1+KmrFDP/KirFDP*(MG/KmrMg)+P/KirP+P/KirP*(MG/KirPMg)+F6P/KirF6P+F6P/KirF6P*(MG/KirF6PMg)+P/KirP*(F6P/KirPF6P)+P/KirP*(F6P/KirPF6P)*(MG/KirPF6PMg)+(FDP-MgFDP)/KirFDP+KdFDPMg/KmrMg*(MgFDP/KirFDP)+AMP/KirAMP+MgFDP/KirFDPMg+MgFDP/KirFDPMg*(MG/KirFDPMgMg)+AMP/KirAMP*((FDP-MgFDP)/KirAMPFDP)))^n)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_514" name="AMP" order="0" role="modifier"/>
        <ParameterDescription key="FunctionParameter_500" name="F6P" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_510" name="FDP" order="2" role="substrate"/>
        <ParameterDescription key="FunctionParameter_512" name="KdFDPMg" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_516" name="KirAMP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_520" name="KirAMPFDP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_502" name="KirF6P" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_504" name="KirF6PMg" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_488" name="KirFDP" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_484" name="KirFDPMg" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_518" name="KirFDPMgMg" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_496" name="KirP" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_506" name="KirPF6P" order="12" role="constant"/>
        <ParameterDescription key="FunctionParameter_508" name="KirPF6PMg" order="13" role="constant"/>
        <ParameterDescription key="FunctionParameter_498" name="KirPMg" order="14" role="constant"/>
        <ParameterDescription key="FunctionParameter_542" name="KitAMP" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_548" name="KitAMPFDP" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_534" name="KitF6P" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_536" name="KitF6PMg" order="18" role="constant"/>
        <ParameterDescription key="FunctionParameter_526" name="KitFDP" order="19" role="constant"/>
        <ParameterDescription key="FunctionParameter_544" name="KitFDPMg" order="20" role="constant"/>
        <ParameterDescription key="FunctionParameter_546" name="KitFDPMgMg" order="21" role="constant"/>
        <ParameterDescription key="FunctionParameter_530" name="KitP" order="22" role="constant"/>
        <ParameterDescription key="FunctionParameter_538" name="KitPF6P" order="23" role="constant"/>
        <ParameterDescription key="FunctionParameter_540" name="KitPF6PMg" order="24" role="constant"/>
        <ParameterDescription key="FunctionParameter_532" name="KitPMg" order="25" role="constant"/>
        <ParameterDescription key="FunctionParameter_486" name="KmrFDP" order="26" role="constant"/>
        <ParameterDescription key="FunctionParameter_492" name="KmrMg" order="27" role="constant"/>
        <ParameterDescription key="FunctionParameter_524" name="KmtFDP" order="28" role="constant"/>
        <ParameterDescription key="FunctionParameter_528" name="KmtMg" order="29" role="constant"/>
        <ParameterDescription key="FunctionParameter_522" name="L0" order="30" role="constant"/>
        <ParameterDescription key="FunctionParameter_490" name="MG" order="31" role="modifier"/>
        <ParameterDescription key="FunctionParameter_463" name="MgFDP" order="32" role="modifier"/>
        <ParameterDescription key="FunctionParameter_494" name="P" order="33" role="product"/>
        <ParameterDescription key="FunctionParameter_483" name="Vmax" order="34" role="constant"/>
        <ParameterDescription key="FunctionParameter_479" name="n" order="35" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_60" name="Function for PPC" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_60">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:43:22Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*n*(PEP*HCO3-OAA*P/Keq)/(KdrPEP*KmrHCO3)/(1+KmrPEP/KdrPEP*(HCO3/KmrHCO3)+KmrOAA/KdrOAA*(P/KmrP)+OAA/KdrOAA+P/KmrP*(OAA/KdrOAA)+HCO3/KmrHCO3*(PEP/KdrPEP)+PEP/KdrPEP)/(1+L0*((1+ACCOA/KeftACCOA+FDP/KeftFDP+FDP/KeftFDP*(ACCOA/KeftFDPACCOA))*(1+KmtPEP/KdtPEP*(HCO3/KmtHCO3)+KmtOAA/KdtOAA*(P/KmtP)+OAA/KdtOAA+P/KmtP*(OAA/KdtOAA)+HCO3/KmtHCO3*(PEP/KdtPEP)+PEP/KdtPEP)*(1+ASP/KeftASP+CYS/KeftCYS+CIT/KeftCIT+FUM/KeftFUM+MAL/KeftMAL+SUC/KeftSUC)/((1+ACCOA/KefrACCOA+FDP/KefrFDP+FDP/KefrFDP*(ACCOA/KefrFDPACCOA))*(1+KmrPEP/KdrPEP*(HCO3/KmrHCO3)+KmrOAA/KdrOAA*(P/KmrP)+OAA/KdrOAA+P/KmrP*(OAA/KdrOAA)+HCO3/KmrHCO3*(PEP/KdrPEP)+PEP/KdrPEP)*(1+ASP/KefrASP+CYS/KefrCYS+CIT/KefrCIT+FUM/KefrFUM+MAL/KefrMAL+SUC/KefrSUC)))^n)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_493" name="ACCOA" order="0" role="modifier"/>
        <ParameterDescription key="FunctionParameter_564" name="ASP" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_572" name="CIT" order="2" role="modifier"/>
        <ParameterDescription key="FunctionParameter_568" name="CYS" order="3" role="modifier"/>
        <ParameterDescription key="FunctionParameter_485" name="FDP" order="4" role="modifier"/>
        <ParameterDescription key="FunctionParameter_576" name="FUM" order="5" role="modifier"/>
        <ParameterDescription key="FunctionParameter_537" name="HCO3" order="6" role="substrate"/>
        <ParameterDescription key="FunctionParameter_505" name="KdrOAA" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_521" name="KdrPEP" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_560" name="KdtOAA" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_554" name="KdtPEP" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_588" name="KefrACCOA" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_594" name="KefrASP" order="12" role="constant"/>
        <ParameterDescription key="FunctionParameter_598" name="KefrCIT" order="13" role="constant"/>
        <ParameterDescription key="FunctionParameter_596" name="KefrCYS" order="14" role="constant"/>
        <ParameterDescription key="FunctionParameter_590" name="KefrFDP" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_592" name="KefrFDPACCOA" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_600" name="KefrFUM" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_602" name="KefrMAL" order="18" role="constant"/>
        <ParameterDescription key="FunctionParameter_604" name="KefrSUC" order="19" role="constant"/>
        <ParameterDescription key="FunctionParameter_489" name="KeftACCOA" order="20" role="constant"/>
        <ParameterDescription key="FunctionParameter_566" name="KeftASP" order="21" role="constant"/>
        <ParameterDescription key="FunctionParameter_574" name="KeftCIT" order="22" role="constant"/>
        <ParameterDescription key="FunctionParameter_570" name="KeftCYS" order="23" role="constant"/>
        <ParameterDescription key="FunctionParameter_469" name="KeftFDP" order="24" role="constant"/>
        <ParameterDescription key="FunctionParameter_550" name="KeftFDPACCOA" order="25" role="constant"/>
        <ParameterDescription key="FunctionParameter_578" name="KeftFUM" order="26" role="constant"/>
        <ParameterDescription key="FunctionParameter_582" name="KeftMAL" order="27" role="constant"/>
        <ParameterDescription key="FunctionParameter_586" name="KeftSUC" order="28" role="constant"/>
        <ParameterDescription key="FunctionParameter_525" name="Keq" order="29" role="constant"/>
        <ParameterDescription key="FunctionParameter_517" name="KmrHCO3" order="30" role="constant"/>
        <ParameterDescription key="FunctionParameter_509" name="KmrOAA" order="31" role="constant"/>
        <ParameterDescription key="FunctionParameter_501" name="KmrP" order="32" role="constant"/>
        <ParameterDescription key="FunctionParameter_513" name="KmrPEP" order="33" role="constant"/>
        <ParameterDescription key="FunctionParameter_556" name="KmtHCO3" order="34" role="constant"/>
        <ParameterDescription key="FunctionParameter_558" name="KmtOAA" order="35" role="constant"/>
        <ParameterDescription key="FunctionParameter_562" name="KmtP" order="36" role="constant"/>
        <ParameterDescription key="FunctionParameter_552" name="KmtPEP" order="37" role="constant"/>
        <ParameterDescription key="FunctionParameter_497" name="L0" order="38" role="constant"/>
        <ParameterDescription key="FunctionParameter_580" name="MAL" order="39" role="modifier"/>
        <ParameterDescription key="FunctionParameter_533" name="OAA" order="40" role="product"/>
        <ParameterDescription key="FunctionParameter_529" name="P" order="41" role="product"/>
        <ParameterDescription key="FunctionParameter_541" name="PEP" order="42" role="substrate"/>
        <ParameterDescription key="FunctionParameter_584" name="SUC" order="43" role="modifier"/>
        <ParameterDescription key="FunctionParameter_549" name="Vmax" order="44" role="constant"/>
        <ParameterDescription key="FunctionParameter_545" name="n" order="45" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_61" name="Function for PCK" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_61">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:54:46Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(MgATP*OAA-HCO3*MgADP*PEP/Keq)/(KmATP*KmOAA)/(1+HCO3/KmHCO3+HCO3/KmHCO3*(ADP/KmADP)+MgADP/KmADP+MgATP/KmATP+OAA/KmOAA+MgATP/KmATP*(OAA/KmOAA)+HCO3/KmHCO3*(PEP/KmPEP)+PEP/KmPEP+HCO3/KmHCO3*(MgADP/KmADP)*(PEP/KmPEP)+MgADP/KmADP*(PEP/KmPEP))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_565" name="ADP" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_593" name="HCO3" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_581" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_561" name="KmADP" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_577" name="KmATP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_569" name="KmHCO3" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_573" name="KmOAA" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_557" name="KmPEP" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_589" name="MgADP" order="8" role="modifier"/>
        <ParameterDescription key="FunctionParameter_601" name="MgATP" order="9" role="modifier"/>
        <ParameterDescription key="FunctionParameter_597" name="OAA" order="10" role="substrate"/>
        <ParameterDescription key="FunctionParameter_585" name="PEP" order="11" role="product"/>
        <ParameterDescription key="FunctionParameter_605" name="Vmax" order="12" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_62" name="Function for PPS" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_62">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:43:41Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(MgATP*PYR-AMP*PEP*P*MG/Keq)/(KmATPMg*KmPYR)/(MgATP/KmATPMg+alpha*(P/KdP)*(MgATP/KmATPMg)+alpha*(AMP/KdAMP)*(MgATP/KmATPMg)+alpha*(P/KdP)*(AMP/KdAMP)*(MgATP/KmATPMg)+alpha*(MG/KdMg)*(P/KmP)*(AMP/KdAMP)*(MgATP/KdATPMgPPS)/(W*(1+MG/KdMg))+MgATP/KmATPMg*(AKG/KefAKG)+(1+MG/KdMg)*(AKG/KefAKG)*(PEP/KmPEP)/W+MgATP/KmATPMg*(OAA/KefOAA)+(1+MG/KdMg)*(OAA/KefOAA)*(PEP/KmPEP)/W+MG/KdMg*(P/KmP)*(AMP/KdAMP)/W+alpha*(P/KdP)*(AMP/KdAMP)*(PEP/KmPEP)/W+alpha*(MG/KdMg)*(P/KmP)*(AMP/KdAMP)*(PEP/KmPEP)/W+alpha*(1+MG/KdMg)*(KmAMP/KdAMP*(P/KmP)*(PEP/KmPEP)+AMP/KdAMP*(PEP/KmPEP))/W+(1+MG/KdMg)*(PYR/KmPYR)+MgATP/KmATPMg*(PYR/KmPYR)+KdADPMg/KdMg*(P/KmP)*(MgADP/KefADP)*(AMP/KdAMP)/(W*(1+MG/KdMg))+(ADP-MgADP)/KefADP*(PYR/KmPYR)+KdATPMg/KdMg*(P/KmP)*(AMP/KdAMP)*(MgATP/KefATP)/(W*(1+MG/KdMg))+(ATP-MgATP)/KefATP*(PYR/KmPYR)+(1+MG/KdMg)*(PEP/KmPEP)/W+alpha*(1+MG/KdMg)*(PEP/KdPEP)*(PYR/KmPYR)+(1+MG/KdMg)*(PYR/KdPYR)*(PEP/KmPEP)/W)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_625" name="ADP" order="0" role="modifier"/>
        <ParameterDescription key="FunctionParameter_607" name="AKG" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_579" name="AMP" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_631" name="ATP" order="3" role="substrate"/>
        <ParameterDescription key="FunctionParameter_619" name="KdADPMg" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_515" name="KdAMP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_627" name="KdATPMg" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_539" name="KdATPMgPPS" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_523" name="KdMg" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_507" name="KdP" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_633" name="KdPEP" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_635" name="KdPYR" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_623" name="KefADP" order="12" role="constant"/>
        <ParameterDescription key="FunctionParameter_609" name="KefAKG" order="13" role="constant"/>
        <ParameterDescription key="FunctionParameter_629" name="KefATP" order="14" role="constant"/>
        <ParameterDescription key="FunctionParameter_615" name="KefOAA" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_551" name="Keq" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_617" name="KmAMP" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_475" name="KmATPMg" order="18" role="constant"/>
        <ParameterDescription key="FunctionParameter_531" name="KmP" order="19" role="constant"/>
        <ParameterDescription key="FunctionParameter_611" name="KmPEP" order="20" role="constant"/>
        <ParameterDescription key="FunctionParameter_491" name="KmPYR" order="21" role="constant"/>
        <ParameterDescription key="FunctionParameter_603" name="MG" order="22" role="modifier"/>
        <ParameterDescription key="FunctionParameter_621" name="MgADP" order="23" role="modifier"/>
        <ParameterDescription key="FunctionParameter_563" name="MgATP" order="24" role="modifier"/>
        <ParameterDescription key="FunctionParameter_613" name="OAA" order="25" role="modifier"/>
        <ParameterDescription key="FunctionParameter_595" name="P" order="26" role="product"/>
        <ParameterDescription key="FunctionParameter_587" name="PEP" order="27" role="product"/>
        <ParameterDescription key="FunctionParameter_571" name="PYR" order="28" role="substrate"/>
        <ParameterDescription key="FunctionParameter_555" name="Vmax" order="29" role="constant"/>
        <ParameterDescription key="FunctionParameter_547" name="W" order="30" role="constant"/>
        <ParameterDescription key="FunctionParameter_499" name="alpha" order="31" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_63" name="Function for MAD" type="UserDefined" reversible="false">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_63">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T23:46:08Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*n*MAL*NAD/(KmrMAL*KirNAD)/(1+KmrNAD/KirNAD*(MAL/KmrMAL)+NAD/KirNAD+MAL/KmrMAL*(NAD/KirNAD))*((MG/KmrMg+MN/KmrMn)/(1+MG/KmrMg+MN/KmrMn))/(1+L0*((1+ASP/KeftASP)*(1+MG/KmtMg+MN/KmtMn)*(1+ATP/KeftATP)*(1+ACCOA/KeftACCOA+COA/KeftCOA)*(1+KmtNAD/KitNAD*(MAL/KmtMAL)+NAD/KitNAD+MAL/KmtMAL*(NAD/KitNAD))/((1+ASP/KefrASP)*(1+MG/KmrMg+MN/KmrMn)*(1+ATP/KefrATP)*(1+ACCOA/KefrACCOA+COA/KefrCOA)*(1+KmrNAD/KirNAD*(MAL/KmrMAL)+NAD/KirNAD+MAL/KmrMAL*(NAD/KirNAD))))^n)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_641" name="ACCOA" order="0" role="modifier"/>
        <ParameterDescription key="FunctionParameter_481" name="ASP" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_637" name="ATP" order="2" role="modifier"/>
        <ParameterDescription key="FunctionParameter_645" name="COA" order="3" role="modifier"/>
        <ParameterDescription key="FunctionParameter_659" name="KefrACCOA" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_655" name="KefrASP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_657" name="KefrATP" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_661" name="KefrCOA" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_643" name="KeftACCOA" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_599" name="KeftASP" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_639" name="KeftATP" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_647" name="KeftCOA" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_616" name="KirNAD" order="12" role="constant"/>
        <ParameterDescription key="FunctionParameter_651" name="KitNAD" order="13" role="constant"/>
        <ParameterDescription key="FunctionParameter_620" name="KmrMAL" order="14" role="constant"/>
        <ParameterDescription key="FunctionParameter_543" name="KmrMg" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_511" name="KmrMn" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_612" name="KmrNAD" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_653" name="KmtMAL" order="18" role="constant"/>
        <ParameterDescription key="FunctionParameter_583" name="KmtMg" order="19" role="constant"/>
        <ParameterDescription key="FunctionParameter_567" name="KmtMn" order="20" role="constant"/>
        <ParameterDescription key="FunctionParameter_649" name="KmtNAD" order="21" role="constant"/>
        <ParameterDescription key="FunctionParameter_495" name="L0" order="22" role="constant"/>
        <ParameterDescription key="FunctionParameter_628" name="MAL" order="23" role="substrate"/>
        <ParameterDescription key="FunctionParameter_608" name="MG" order="24" role="modifier"/>
        <ParameterDescription key="FunctionParameter_527" name="MN" order="25" role="modifier"/>
        <ParameterDescription key="FunctionParameter_624" name="NAD" order="26" role="substrate"/>
        <ParameterDescription key="FunctionParameter_636" name="Vmax" order="27" role="constant"/>
        <ParameterDescription key="FunctionParameter_632" name="n" order="28" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_64" name="Function for PDH" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_64">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T23:51:24Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(COA*NAD*PYR-ACCOA*NADH*HCO3/Keq)/(KmCOA*KmNAD*KmPYR)/(ACCOA/KmACCOA+NADH/KmNADH+ACCOA/KmACCOA*(NADH/KmNADH)+COA/KmCOA*(NADH/KmNADH)+ACCOA/KmACCOA*(COA/KmCOA)*(NADH/KmNADH)+NAD/KmNAD*(NADH/KmNADH)+COA/KmCOA*(NAD/KmNAD)*(NADH/KmNADH)+ACCOA/KmACCOA*(PYR/KmPYR)+ACCOA/KmACCOA*(COA/KmCOA)*(PYR/KmPYR)+COA/KmCOA*(1+NAD/KmNAD)*(PYR/KmPYR)+NAD/KmNAD*(1+COA/KmCOA+PYR/KmPYR))/(1+HCO3/KmHCO3)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_646" name="ACCOA" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_658" name="COA" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_638" name="HCO3" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_575" name="Keq" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_610" name="KmACCOA" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_553" name="KmCOA" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_626" name="KmHCO3" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_503" name="KmNAD" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_618" name="KmNADH" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_535" name="KmPYR" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_654" name="NAD" order="10" role="substrate"/>
        <ParameterDescription key="FunctionParameter_642" name="NADH" order="11" role="product"/>
        <ParameterDescription key="FunctionParameter_650" name="PYR" order="12" role="substrate"/>
        <ParameterDescription key="FunctionParameter_662" name="Vmax" order="13" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_65" name="Function for GLT" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_65">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T14:30:02Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(ACCOA*OAA-CIT*COA/Keq)/(KdACCOA0*KmOAA0)/(ACCOA/KdACCOA0*(1+AKG/Ki1AKG+NADH/Ki1NADH)+ACCOA/KdACCOA0*(OAA/KmOAA0)*(1+AKG/Ki2AKG+NADH/Ki2NADH)+(1+ATP/KiATP)*(1+KmACCOA0/KdACCOA0*(OAA/KmOAA0))+KmcsCOA/KdcsCOA*(CIT/KmcsCIT)+KmcsCOA/KdcsCOA*(ACCOA/KdACCOA0)*(CIT/KmcsCIT)+COA/KdcsCOA+CIT/KmcsCIT*(OAA/KmOAA0)+ACCOA/KdACCOA0*(CIT/KdcsCIT)*(OAA/KmOAA0)+KmACCOA0/KdACCOA0*(COA/KdcsCOA)*(OAA/KmOAA0)+CIT/KmcsCIT*(COA/KdcsCOA)*(OAA/KdcsOAA))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_614" name="ACCOA" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_664" name="AKG" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_676" name="ATP" order="2" role="modifier"/>
        <ParameterDescription key="FunctionParameter_591" name="CIT" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_640" name="COA" order="4" role="product"/>
        <ParameterDescription key="FunctionParameter_656" name="KdACCOA0" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_688" name="KdcsCIT" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_684" name="KdcsCOA" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_690" name="KdcsOAA" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_648" name="Keq" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_666" name="Ki1AKG" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_670" name="Ki1NADH" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_672" name="Ki2AKG" order="12" role="constant"/>
        <ParameterDescription key="FunctionParameter_674" name="Ki2NADH" order="13" role="constant"/>
        <ParameterDescription key="FunctionParameter_678" name="KiATP" order="14" role="constant"/>
        <ParameterDescription key="FunctionParameter_680" name="KmACCOA0" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_634" name="KmOAA0" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_686" name="KmcsCIT" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_682" name="KmcsCOA" order="18" role="constant"/>
        <ParameterDescription key="FunctionParameter_668" name="NADH" order="19" role="modifier"/>
        <ParameterDescription key="FunctionParameter_519" name="OAA" order="20" role="substrate"/>
        <ParameterDescription key="FunctionParameter_630" name="Vmax" order="21" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_66" name="Function for ACN_1" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_66">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-13T10:11:09Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(CIT-ACO/Keq)/KmCIT/(1+CIT/KmCIT+ACO/KmACO+ICIT/KmICIT)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_683" name="ACO" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_687" name="CIT" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_667" name="ICIT" order="2" role="modifier"/>
        <ParameterDescription key="FunctionParameter_679" name="Keq" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_671" name="KmACO" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_675" name="KmCIT" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_663" name="KmICIT" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_691" name="Vmax" order="7" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_67" name="Function for ACN_2" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_67">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-13T10:11:14Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(ACO-ICIT/Keq)/KmACO/(1+ACO/KmACO+ICIT/KmICIT+CIT/KmCIT)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_669" name="ACO" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_606" name="CIT" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_677" name="ICIT" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_685" name="Keq" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_652" name="KmACO" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_692" name="KmCIT" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_559" name="KmICIT" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_660" name="Vmax" order="7" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_68" name="Function for ICD" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_68">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T09:24:11Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        icd*kcat*(ICIT*NADP-AKG*NADPH/Keq)/(KmICIT*KmNADP)/((1+ICIT/KmICIT)*(1+NADP/KmNADP)+(1+AKG/KmAKG)*(1+NADPH/KmNADPH)-1)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_694" name="AKG" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_689" name="ICIT" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_698" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_704" name="KmAKG" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_700" name="KmICIT" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_702" name="KmNADP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_706" name="KmNADPH" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_673" name="NADP" order="7" role="substrate"/>
        <ParameterDescription key="FunctionParameter_696" name="NADPH" order="8" role="product"/>
        <ParameterDescription key="FunctionParameter_693" name="icd" order="9" role="modifier"/>
        <ParameterDescription key="FunctionParameter_487" name="kcat" order="10" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_69" name="Function for LPD" type="UserDefined" reversible="false">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_69">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T23:52:31Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*COA*AKG*NAD*(1-AKG/KdAKG)/(KmCOA*KmAKG*KmNAD)/(COA/KmCOA*(AKG/KmAKG)+COA/KmCOA*(NAD/KmNAD)+AKG/KmAKG*(NAD/KmNAD)+COA/KmCOA*(AKG/KmAKG)*(NAD/KmNAD)-AKG/KdAKG*(COA/KmCOA*(AKG/KmAKG)+AKG/KmAKG*(NAD/KmNAD)+alpha*(COA/KmCOA)*(AKG/KmAKG)*(NAD/KmNAD)))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_699" name="AKG" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_703" name="COA" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_681" name="KdAKG" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_709" name="KmAKG" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_622" name="KmCOA" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_711" name="KmNAD" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_695" name="NAD" order="6" role="substrate"/>
        <ParameterDescription key="FunctionParameter_707" name="Vmax" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_713" name="alpha" order="8" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_70" name="Function for SK" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_70">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:44:49Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(ADP*SUCCOA*P-ATP*COA*SUC/Keq)/(KmADP*KmSUCCOA*KmP)/((1+ADP/KmADP)*(1+SUCCOA/KmSUCCOA)*(1+P/KmP)+(1+ATP/KmATP)*(1+COA/KmCOA)*(1+SUC/KmSUC)-1)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_710" name="ADP" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_705" name="ATP" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_716" name="COA" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_720" name="Keq" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_722" name="KmADP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_728" name="KmATP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_730" name="KmCOA" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_726" name="KmP" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_732" name="KmSUC" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_724" name="KmSUCCOA" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_697" name="P" order="10" role="substrate"/>
        <ParameterDescription key="FunctionParameter_718" name="SUC" order="11" role="product"/>
        <ParameterDescription key="FunctionParameter_644" name="SUCCOA" order="12" role="substrate"/>
        <ParameterDescription key="FunctionParameter_714" name="Vmax" order="13" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_71" name="Function for SDH" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_71">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T13:52:00Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(SUC*Q-FUM*QH2/Keq)/(KefSUC*KmQ)/(1+FUM/KefFUM+KmSUC/KefSUC*(Q/KmQ)+KmFUM/KefFUM*(QH2/KmQH2)+FUM/KefFUM*(QH2/KmQH2)+SUC/KefSUC+SUC/KefSUC*(Q/KmQ))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_721" name="FUM" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_736" name="KefFUM" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_708" name="KefSUC" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_701" name="Keq" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_740" name="KmFUM" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_734" name="KmQ" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_742" name="KmQH2" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_738" name="KmSUC" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_725" name="Q" order="8" role="substrate"/>
        <ParameterDescription key="FunctionParameter_717" name="QH2" order="9" role="product"/>
        <ParameterDescription key="FunctionParameter_729" name="SUC" order="10" role="substrate"/>
        <ParameterDescription key="FunctionParameter_733" name="Vmax" order="11" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_72" name="Function for FUMA" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_72">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-11T17:11:24Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(FUM-MAL/Keq)/KmFUM/(1+FUM/KmFUM+MAL/KmMAL)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_739" name="FUM" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_665" name="Keq" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_719" name="KmFUM" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_727" name="KmMAL" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_735" name="MAL" order="4" role="product"/>
        <ParameterDescription key="FunctionParameter_743" name="Vmax" order="5" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_73" name="Function for MQO" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_73">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T17:17:58Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(MAL*Q-OAA*QH2/Keq)/(KmMAL*KmQ)/((1+MAL/KmMAL)*(1+Q/KmQ)+(1+OAA/KmOAA)*(1+QH2/KmQH2)-1)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_748" name="Keq" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_750" name="KmMAL" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_754" name="KmOAA" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_752" name="KmQ" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_756" name="KmQH2" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_715" name="MAL" order="5" role="substrate"/>
        <ParameterDescription key="FunctionParameter_744" name="OAA" order="6" role="product"/>
        <ParameterDescription key="FunctionParameter_737" name="Q" order="7" role="substrate"/>
        <ParameterDescription key="FunctionParameter_746" name="QH2" order="8" role="product"/>
        <ParameterDescription key="FunctionParameter_731" name="Vmax" order="9" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_74" name="Function for MDH" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_74">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T23:52:49Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(NADH*OAA-MAL*NAD/Keq)/(KiNADH*KmOAA)/(1+KmNAD/KiNAD*(MAL/KmMAL)+NAD/KiNAD+MAL/KmMAL*(NAD/KiNAD)+NADH/KiNADH+KmNAD/KiNAD*(MAL/KmMAL)*(NADH/KiNADH)+KmNADH/KiNADH*(OAA/KmOAA)+KmNADH/KiNADH*(NAD/KiNAD)*(OAA/KmOAA)+MAL*NAD*OAA/(KiNAD*KiOAA*KmMAL)+NADH/KiNADH*(OAA/KmOAA)+Keq*(KiNADH*KmOAA/(KiNAD*KmMAL))*(KmNAD/KiNAD)*(MAL/KmMAL)*(NADH/KmNADH)*(OAA/KiOAA))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_758" name="Keq" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_766" name="KiNAD" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_760" name="KiNADH" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_772" name="KiOAA" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_768" name="KmMAL" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_764" name="KmNAD" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_770" name="KmNADH" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_762" name="KmOAA" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_745" name="MAL" order="8" role="product"/>
        <ParameterDescription key="FunctionParameter_712" name="NAD" order="9" role="modifier"/>
        <ParameterDescription key="FunctionParameter_753" name="NADH" order="10" role="modifier"/>
        <ParameterDescription key="FunctionParameter_749" name="OAA" order="11" role="substrate"/>
        <ParameterDescription key="FunctionParameter_757" name="Vmax" order="12" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_75" name="Function for ACEA" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_75">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T10:00:43Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(ICIT-GLX*SUC/Keq)/KmICIT/(1+ICIT/KmICIT*(1+PEP/KdPEPicit)+SUC/KdSUC*(1+ICIT/KdICITsuc)+KmSUC/KdSUC*(GLX/KmGLX)*(1+PEP/KdPEPglx)+GLX/KmGLX*(SUC/KdSUC)+PEP/KdPEP+PGA3/KdPGA3)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_765" name="GLX" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_769" name="ICIT" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_779" name="KdICITsuc" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_787" name="KdPEP" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_785" name="KdPEPglx" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_775" name="KdPEPicit" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_791" name="KdPGA3" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_777" name="KdSUC" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_723" name="Keq" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_783" name="KmGLX" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_747" name="KmICIT" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_781" name="KmSUC" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_755" name="PEP" order="12" role="modifier"/>
        <ParameterDescription key="FunctionParameter_789" name="PGA3" order="13" role="modifier"/>
        <ParameterDescription key="FunctionParameter_761" name="SUC" order="14" role="product"/>
        <ParameterDescription key="FunctionParameter_773" name="Vmax" order="15" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_76" name="Function for ACEB" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_76">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T14:32:44Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(ACCOA*GLX-COA*MAL/Keq)/(KmACCOA*KmGLX)/((1+ACCOA/KmACCOA)*(1+GLX/KmGLX)+(1+COA/KmCOA)*(1+MAL/KmMAL)-1)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_788" name="ACCOA" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_780" name="COA" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_784" name="GLX" order="2" role="substrate"/>
        <ParameterDescription key="FunctionParameter_751" name="Keq" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_759" name="KmACCOA" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_793" name="KmCOA" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_767" name="KmGLX" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_795" name="KmMAL" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_776" name="MAL" order="8" role="product"/>
        <ParameterDescription key="FunctionParameter_792" name="Vmax" order="9" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_77" name="Function for ACEK_1" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_77">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T13:25:46Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        k*(ATP*icd-ADP*icdP/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_778" name="ADP" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_771" name="ATP" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_797" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_741" name="icd" order="3" role="substrate"/>
        <ParameterDescription key="FunctionParameter_786" name="icdP" order="4" role="product"/>
        <ParameterDescription key="FunctionParameter_796" name="k" order="5" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_78" name="Function for ACEK_2" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_78">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T14:20:35Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        k*(icdP-icd*P/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_801" name="Keq" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_799" name="P" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_763" name="icd" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_782" name="icdP" order="3" role="substrate"/>
        <ParameterDescription key="FunctionParameter_798" name="k" order="4" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_79" name="Function for EDD" type="UserDefined" reversible="true">
      <Expression>
        Vmax*(PGN-KDPG/Keq)/KmPGN/(1+PGN/KmPGN+KDPG/KmKDPG)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_790" name="KDPG" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_804" name="Keq" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_808" name="KmKDPG" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_806" name="KmPGN" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_794" name="PGN" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_802" name="Vmax" order="5" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_80" name="Function for EDA" type="UserDefined" reversible="true">
      <Expression>
        Vmax*(KDPG-GAP*PYR/Keq)/KmKDPG/(1+KDPG/KmKDPG+(1+GAP/KmGAP)*(1+PYR/KmPYR)-1)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_774" name="GAP" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_805" name="KDPG" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_812" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_816" name="KmGAP" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_814" name="KmKDPG" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_818" name="KmPYR" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_810" name="PYR" order="6" role="product"/>
        <ParameterDescription key="FunctionParameter_809" name="Vmax" order="7" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_81" name="Function for PNT" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_81">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T13:50:16Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        k*(NAD*NADPH-NADH*NADP/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_822" name="Keq" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_815" name="NAD" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_803" name="NADH" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_820" name="NADP" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_811" name="NADPH" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_819" name="k" order="5" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_82" name="Function for ADK" type="UserDefined" reversible="true">
      <Expression>
        k*(AMP*ATP-ADP^2/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_824" name="ADP" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_807" name="AMP" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_813" name="ATP" order="2" role="substrate"/>
        <ParameterDescription key="FunctionParameter_826" name="Keq" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_823" name="k" order="4" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_83" name="Function for CYA" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_83">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-27T11:43:11Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        k*(ATP-CAMP*P^2/Keq)*eiiaP/(eiiaP+KaeiiaP)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_817" name="ATP" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_821" name="CAMP" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_831" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_829" name="P" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_827" name="k" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_833" name="eiiaP" order="5" role="modifier"/>
        <ParameterDescription key="FunctionParameter_835" name="KaeiiaP" order="6" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_84" name="Function for DOS" type="UserDefined" reversible="true">
      <Expression>
        k*(CAMP-AMP/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_828" name="AMP" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_832" name="CAMP" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_825" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_836" name="k" order="3" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_85" name="Function for ACK" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_85">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:42:30Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        cell*Vmax*(ACP*ADP-ACEx*ATP/Keq)/(KmACP*KmADP)/((1+ACP/KmACP+ACEx/KmACE)*(1+ADP/KmADP+ATP/KmATP))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_842" name="ACEx" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_838" name="ACP" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_840" name="ADP" order="2" role="substrate"/>
        <ParameterDescription key="FunctionParameter_844" name="ATP" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_846" name="Keq" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_852" name="KmACE" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_848" name="KmACP" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_850" name="KmADP" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_854" name="KmATP" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_830" name="Vmax" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_837" name="cell" order="10" role="volume"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_86" name="Function for ACS" type="UserDefined" reversible="false">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_86">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T14:33:58Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        cell*Vmax*ACEx*ATP*COA/(KmACE*KmATP*KmCOA)/((1+ACEx/KmACE)*(1+ATP/KmATP)*(1+COA/KmCOA))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_847" name="ACEx" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_843" name="ATP" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_839" name="COA" order="2" role="substrate"/>
        <ParameterDescription key="FunctionParameter_800" name="KmACE" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_857" name="KmATP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_859" name="KmCOA" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_851" name="Vmax" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_855" name="cell" order="7" role="volume"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_87" name="Function for PTA" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_87">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:44:15Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(ACCOA*P-ACP*COA/Keq)/(KiACCOA*KmP)/(1+ACCOA/KiACCOA+P/KiP+ACP/KiACP+COA/KiCOA+ACCOA*P/(KiACCOA*KmP)+ACP*COA/(KmACP*KiCOA))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_856" name="ACCOA" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_849" name="ACP" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_861" name="COA" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_863" name="Keq" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_865" name="KiACCOA" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_871" name="KiACP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_873" name="KiCOA" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_869" name="KiP" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_875" name="KmACP" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_867" name="KmP" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_841" name="P" order="10" role="substrate"/>
        <ParameterDescription key="FunctionParameter_860" name="Vmax" order="11" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_88" name="Function for PTS_0" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_88">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-22T14:03:00Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        kF*ei*PEP^2/(KmPEP^2+PEP^2)-kR*eiP*PYR^2/(KmPYR^2+PYR^2)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_864" name="KmPEP" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_879" name="KmPYR" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_868" name="PEP" order="2" role="substrate"/>
        <ParameterDescription key="FunctionParameter_877" name="PYR" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_872" name="ei" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_834" name="eiP" order="5" role="product"/>
        <ParameterDescription key="FunctionParameter_876" name="kF" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_853" name="kR" order="7" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_89" name="Function for PTS_4" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_89">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-22T14:03:01Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        cell*(kF*eiicbP*GLCx/(KmGLC+GLCx)-kR*eiicb*G6P/(KmG6P+G6P))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_887" name="G6P" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_870" name="GLCx" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_889" name="KmG6P" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_881" name="KmGLC" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_880" name="cell" order="4" role="volume"/>
        <ParameterDescription key="FunctionParameter_885" name="eiicb" order="5" role="product"/>
        <ParameterDescription key="FunctionParameter_862" name="eiicbP" order="6" role="substrate"/>
        <ParameterDescription key="FunctionParameter_858" name="kF" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_883" name="kR" order="8" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_90" name="Function for GLC_FEED" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_90">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2014-12-29T21:58:51Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        cell*FEED/extracellular
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_886" name="FEED" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_890" name="cell" order="1" role="volume"/>
        <ParameterDescription key="FunctionParameter_882" name="extracellular" order="2" role="volume"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_91" name="Function for ATP_MAINTENANCE" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_91">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-14T11:56:01Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(ATP-ADP*P/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_874" name="Vmax" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_845" name="ADP" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_891" name="P" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_888" name="ATP" order="3" role="substrate"/>
        <ParameterDescription key="FunctionParameter_893" name="Keq" order="4" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_92" name="Function for ATP_SYN" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_92">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T14:25:01Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*((log(Hout/Hin)/log(10))^4)/(1+((log(Hout/Hin)/log(10))^4))*(ADP*P-ATP/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_894" name="Vmax" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_884" name="Hin" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_878" name="Hout" order="2" role="modifier"/>
        <ParameterDescription key="FunctionParameter_896" name="ADP" order="3" role="substrate"/>
        <ParameterDescription key="FunctionParameter_898" name="P" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_900" name="ATP" order="5" role="product"/>
        <ParameterDescription key="FunctionParameter_902" name="Keq" order="6" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_93" name="Function for NDH1" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_93">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T14:28:44Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax/(1+((log(Hout/Hin)/log(10))^2))*(NADH*Q-NAD*QH2/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_903" name="Vmax" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_895" name="Hin" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_899" name="Hout" order="2" role="modifier"/>
        <ParameterDescription key="FunctionParameter_892" name="NADH" order="3" role="substrate"/>
        <ParameterDescription key="FunctionParameter_905" name="Q" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_907" name="NAD" order="5" role="product"/>
        <ParameterDescription key="FunctionParameter_909" name="QH2" order="6" role="product"/>
        <ParameterDescription key="FunctionParameter_911" name="Keq" order="7" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_94" name="Function for CYTBO" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_94">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T14:29:37Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax/(1+((log(Hout/Hin)/log(10))^2))*(QH2^2*O2-Q^2/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_912" name="Vmax" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_904" name="Hin" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_908" name="Hout" order="2" role="modifier"/>
        <ParameterDescription key="FunctionParameter_917" name="Keq" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_913" name="O2" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_897" name="QH2" order="5" role="substrate"/>
        <ParameterDescription key="FunctionParameter_915" name="Q" order="6" role="product"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_95" name="Function for SQR" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_95">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T14:33:09Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(FADH2*Q-FAD*QH2/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_914" name="FADH2" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_866" name="Q" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_910" name="FAD" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_920" name="QH2" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_922" name="Keq" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_918" name="Vmax" order="5" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_96" name="Function for NDH2" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_96">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-17T09:50:33Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(NADH*Q-NAD*QH2/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_923" name="Vmax" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_919" name="NADH" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_901" name="Q" order="2" role="substrate"/>
        <ParameterDescription key="FunctionParameter_924" name="NAD" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_926" name="QH2" order="4" role="product"/>
        <ParameterDescription key="FunctionParameter_928" name="Keq" order="5" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_97" name="Function for GROWTH" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_97">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-21T14:35:31Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax * G6P * E4P * PGA3 * OAA * AKG * PYR * R5P * PEP * GAP * F6P * NADPH * ACCOA * NAD * ATP / (KmG6P * KmE4P * KmPGA3 * KmOAA * KmAKG * KmPYR * KmR5P * KmPEP * KmGAP * KmF6P * KmNADPH * KmACCOA * KmNAD * KmATP) / ((1+G6P/KmG6P)*(1+E4P/KmE4P)*(1+PGA3/KmPGA3)*(1+OAA/KmOAA)*(1+AKG/KmAKG)*(1+PYR/KmPYR)*(1+R5P/KmR5P)*(1+PEP/KmPEP)*(1+GAP/KmGAP)*(1+F6P/KmF6P)*(1+NADPH/KmNADPH)*(1+ACCOA/KmACCOA)*(1+NAD/KmNAD)*(1+ATP/KmATP))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_929" name="Vmax" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_925" name="G6P" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_906" name="E4P" order="2" role="substrate"/>
        <ParameterDescription key="FunctionParameter_930" name="PGA3" order="3" role="substrate"/>
        <ParameterDescription key="FunctionParameter_932" name="OAA" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_934" name="AKG" order="5" role="substrate"/>
        <ParameterDescription key="FunctionParameter_936" name="PYR" order="6" role="substrate"/>
        <ParameterDescription key="FunctionParameter_938" name="R5P" order="7" role="substrate"/>
        <ParameterDescription key="FunctionParameter_940" name="PEP" order="8" role="substrate"/>
        <ParameterDescription key="FunctionParameter_942" name="GAP" order="9" role="substrate"/>
        <ParameterDescription key="FunctionParameter_944" name="F6P" order="10" role="substrate"/>
        <ParameterDescription key="FunctionParameter_946" name="NADPH" order="11" role="substrate"/>
        <ParameterDescription key="FunctionParameter_948" name="ACCOA" order="12" role="substrate"/>
        <ParameterDescription key="FunctionParameter_950" name="NAD" order="13" role="substrate"/>
        <ParameterDescription key="FunctionParameter_952" name="ATP" order="14" role="substrate"/>
        <ParameterDescription key="FunctionParameter_954" name="KmG6P" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_956" name="KmE4P" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_958" name="KmPGA3" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_960" name="KmOAA" order="18" role="constant"/>
        <ParameterDescription key="FunctionParameter_962" name="KmAKG" order="19" role="constant"/>
        <ParameterDescription key="FunctionParameter_964" name="KmPYR" order="20" role="constant"/>
        <ParameterDescription key="FunctionParameter_966" name="KmR5P" order="21" role="constant"/>
        <ParameterDescription key="FunctionParameter_968" name="KmPEP" order="22" role="constant"/>
        <ParameterDescription key="FunctionParameter_970" name="KmGAP" order="23" role="constant"/>
        <ParameterDescription key="FunctionParameter_972" name="KmF6P" order="24" role="constant"/>
        <ParameterDescription key="FunctionParameter_974" name="KmNADPH" order="25" role="constant"/>
        <ParameterDescription key="FunctionParameter_976" name="KmACCOA" order="26" role="constant"/>
        <ParameterDescription key="FunctionParameter_978" name="KmNAD" order="27" role="constant"/>
        <ParameterDescription key="FunctionParameter_980" name="KmATP" order="28" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_98" name="Function for PIT" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_98">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-31T12:39:10Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        (Vmax*(((log(Hout/Hin)/log(10))^2/(1+(log(Hout/Hin)/log(10))^2))*(Pp/(KmPp+Pp))-(Kr/(1+(log(Hout/Hin)/log(10))^2))*(P/(KmP+P))))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_981" name="Vmax" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_977" name="Hout" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_973" name="Hin" order="2" role="modifier"/>
        <ParameterDescription key="FunctionParameter_969" name="Pp" order="3" role="substrate"/>
        <ParameterDescription key="FunctionParameter_965" name="KmPp" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_961" name="Kr" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_957" name="P" order="6" role="product"/>
        <ParameterDescription key="FunctionParameter_953" name="KmP" order="7" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_99" name="Function for GL6P_HYDROLYSIS" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_99">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-22T11:04:41Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        KGl6Phydrol*(GL6P-PGN/KeqGl6Phydrol)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_951" name="KGl6Phydrol" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_959" name="GL6P" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_967" name="PGN" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_975" name="KeqGl6Phydrol" order="3" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_100" name="Function for XCH_RMM" type="UserDefined" reversible="unspecified">
      <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_100">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2016-05-27T10:08:32Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
      </MiriamAnnotation>
      <Expression>
        Vmax*(S/Km-P/Km)/(1+S/Km+P/Km)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_979" name="Vmax" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_963" name="S" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_945" name="P" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_949" name="Km" order="3" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
  </ListOfFunctions>
  <Model key="Model_3" name="Millard2016 - E. coli central carbon and energy metabolism" simulationType="time" timeUnit="s" volumeUnit="l" areaUnit="m²" lengthUnit="m" quantityUnit="mmol" type="deterministic" avogadroConstant="6.02214179e+023">
    <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#">
  <rdf:Description rdf:about="#Model_3">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-05-10T22:00:37Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
    <dcterms:creator>
      <rdf:Description>
        <vCard:EMAIL>millard@insa-toulouse.fr</vCard:EMAIL>
        <vCard:N>
          <rdf:Description>
            <vCard:Family>Millard</vCard:Family>
            <vCard:Given>Pierre</vCard:Given>
          </rdf:Description>
        </vCard:N>
        <vCard:ORG>
          <rdf:Description>
            <vCard:Orgname>INRA</vCard:Orgname>
          </rdf:Description>
        </vCard:ORG>
      </rdf:Description>
    </dcterms:creator>
    <dcterms:creator>
      <rdf:Description>
        <vCard:EMAIL>kieran.smallbone@manchester.ac.uk</vCard:EMAIL>
        <vCard:N>
          <rdf:Description>
            <vCard:Family>Smallbone</vCard:Family>
            <vCard:Given>Kieran</vCard:Given>
          </rdf:Description>
        </vCard:N>
        <vCard:ORG>
          <rdf:Description>
            <vCard:Orgname>University of Manchester</vCard:Orgname>
          </rdf:Description>
        </vCard:ORG>
      </rdf:Description>
    </dcterms:creator>
    <dcterms:creator>
      <rdf:Description>
        <vCard:EMAIL>pedro.mendes@manchester.ac.uk</vCard:EMAIL>
        <vCard:N>
          <rdf:Description>
            <vCard:Family>Mendes</vCard:Family>
            <vCard:Given>Pedro</vCard:Given>
          </rdf:Description>
        </vCard:N>
        <vCard:ORG>
          <rdf:Description>
            <vCard:Orgname>University of Manchester</vCard:Orgname>
          </rdf:Description>
        </vCard:ORG>
      </rdf:Description>
    </dcterms:creator>
  </rdf:Description>
</rdf:RDF>

    </MiriamAnnotation>
    <Comment>
      This model represents the central carbon and energy metabolic network of Escherichia coli growing under glucose limitation.
    </Comment>
    <ListOfCompartments>
      <Compartment key="Compartment_1" name="cell_cytoplasm" simulationType="fixed" dimensionality="3">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Compartment_1">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T16:22:10Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Compartment>
      <Compartment key="Compartment_3" name="extracellular" simulationType="fixed" dimensionality="3">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Compartment_3">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T16:22:10Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Compartment>
      <Compartment key="Compartment_5" name="cell_periplasm" simulationType="fixed" dimensionality="3">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Compartment_5">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-16T10:08:42Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Compartment>
    </ListOfCompartments>
    <ListOfMetabolites>
      <Metabolite key="Metabolite_1" name="ACCOA" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_1">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-17T08:09:06Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_3" name="ACO" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_3">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-14T12:22:27Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_5" name="ACP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_5">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-14T12:22:29Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_7" name="AKG" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_7">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-17T08:08:54Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_9" name="BPG" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_9">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T14:57:13Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_11" name="CIT" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" />
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_13" name="DAP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_13">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T10:03:10Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_15" name="E4P" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_15">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-16T11:41:38Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_17" name="F6P" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_17">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T16:25:24Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_19" name="FDP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_19">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:49:14Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_21" name="FUM" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" />
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_23" name="G6P" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" />
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_25" name="GAP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_25">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-03T00:12:37Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_27" name="GL6P" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_27">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:53:03Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_29" name="GLX" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_29">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:53:08Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_31" name="ICIT" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_31">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T14:59:12Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_33" name="KDPG" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_33">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-24T10:12:51Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_35" name="MAL" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" />
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_37" name="NAD" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_37">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T23:44:07Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_39" name="NADH" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_39">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T23:44:56Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_41" name="NADP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_41">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T17:09:49Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_43" name="NADPH" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_43">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:58:22Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_45" name="OAA" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_45">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:45:19Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_47" name="PEP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_47">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-24T16:37:45Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_49" name="PGA2" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" />
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_51" name="PGA3" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" />
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_53" name="PGN" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" />
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_55" name="PYR" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_55">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-31T12:41:19Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_57" name="Q" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_57">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-14T12:12:52Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_59" name="QH2" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_59">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-14T12:12:53Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_61" name="R5P" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" />
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_63" name="RU5P" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" />
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_65" name="S7P" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_65">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-24T10:25:24Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_67" name="SUC" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_67">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:46:56Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_69" name="SUCCOA" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_69">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-14T13:47:22Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_71" name="X5P" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_71">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-22T20:51:59Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_73" name="ei" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_73">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T14:56:48Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_75" name="eiP" simulationType="reactions" compartment="Compartment_1">
      </Metabolite>
      <Metabolite key="Metabolite_77" name="eiia" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_77">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-16T09:44:24Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_79" name="eiiaP" simulationType="reactions" compartment="Compartment_1">
      </Metabolite>
      <Metabolite key="Metabolite_81" name="eiicb" simulationType="reactions" compartment="Compartment_1">
      </Metabolite>
      <Metabolite key="Metabolite_83" name="eiicbP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_83">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T11:48:22Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_85" name="hpr" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_85">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:47:23Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_87" name="hprP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_87">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-10T10:15:23Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_89" name="icd" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_89">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:53:36Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_91" name="icdP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_91">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:53:40Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_93" name="tal" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_93">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T10:37:13Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_95" name="talC3" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_95">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T10:37:13Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_97" name="tkt" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_97">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T10:37:48Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_99" name="tktC2" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_99">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T10:37:48Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_101" name="ADP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_101">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:47:02Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_103" name="AMP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_103">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:47:05Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_105" name="ATP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_105">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:47:08Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_107" name="CAMP" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_107">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:47:11Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_109" name="COA" simulationType="fixed" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_109">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T14:35:06Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_111" name="HCO3" simulationType="fixed" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_111">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:47:23Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_113" name="P" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_113">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:45:13Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_115" name="MG" simulationType="fixed" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_115">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:50:03Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_117" name="MgADP" simulationType="assignment" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_117">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T00:28:29Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <Expression>
          &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG],Reference=Concentration&gt;*&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ADP],Reference=Concentration&gt;/(&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdADPMg],Reference=Value&gt;+&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG],Reference=Concentration&gt;)
        </Expression>
      </Metabolite>
      <Metabolite key="Metabolite_119" name="MgATP" simulationType="assignment" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_119">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-23T11:39:48Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <Expression>
          &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG],Reference=Concentration&gt;*&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ATP],Reference=Concentration&gt;/(&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdATPMg],Reference=Value&gt;+&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG],Reference=Concentration&gt;)
        </Expression>
      </Metabolite>
      <Metabolite key="Metabolite_121" name="MgFDP" simulationType="assignment" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_121">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T00:28:01Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <Expression>
          &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG],Reference=Concentration&gt;*&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[FDP],Reference=Concentration&gt;/(&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdFDPMg],Reference=Value&gt;+&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG],Reference=Concentration&gt;)
        </Expression>
      </Metabolite>
      <Metabolite key="Metabolite_123" name="ASP" simulationType="fixed" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_123">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-14T12:17:44Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_125" name="CYS" simulationType="fixed" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_125">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-16T14:03:20Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_127" name="MN" simulationType="fixed" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_127">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-27T11:44:38Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_129" name="Hin" simulationType="fixed" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_129">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T10:15:49Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_131" name="H2O" simulationType="fixed" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_131">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T10:15:34Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_133" name="O2" simulationType="fixed" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_133">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T10:59:46Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_135" name="FAD" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_135">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T15:23:15Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_137" name="FADH2" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_137">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T15:23:13Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_139" name="ACE" simulationType="reactions" compartment="Compartment_1">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_139">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-14T12:22:26Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_141" name="GLCx" simulationType="reactions" compartment="Compartment_3">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_141">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:49:09Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_143" name="Px" simulationType="fixed" compartment="Compartment_3">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_143">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-24T10:32:17Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_145" name="ACEx" simulationType="reactions" compartment="Compartment_3">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_145">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-03T00:04:02Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_147" name="Hout" simulationType="reactions" compartment="Compartment_5">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_147">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T10:15:59Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_149" name="GLCp" simulationType="reactions" compartment="Compartment_5">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_149">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-23T16:11:40Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_151" name="Pp" simulationType="reactions" compartment="Compartment_5">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_151">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-31T12:41:29Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_153" name="ACEp" simulationType="reactions" compartment="Compartment_5">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_153">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-03T00:12:13Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </Metabolite>
    </ListOfMetabolites>
    <ListOfModelValues>
      <ModelValue key="ModelValue_0" name="FEED" simulationType="fixed">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_0">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T10:22:21Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_1" name="KdADPMg" simulationType="fixed">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_1">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T19:54:50Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_2" name="KdATPMg" simulationType="fixed">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_2">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-20T18:47:20Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_3" name="KdFDPMg" simulationType="fixed">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_3">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-20T18:47:20Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_4" name="KmICIT_ACN" simulationType="fixed">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_4">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-13T02:18:45Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_5" name="KmCIT_ACN" simulationType="fixed">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_5">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-13T02:18:59Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_6" name="KmACO_ACN" simulationType="fixed">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_6">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-13T02:19:05Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_7" name="KeqNDH" simulationType="fixed">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_7">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-25T10:29:33Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
    </ListOfModelValues>
    <ListOfReactions>
      <Reaction key="Reaction_0" name="PGI" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_0">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:08:22Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_23" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_17" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_47" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_53" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4344" name="Keq" value="0.36"/>
          <Constant key="Parameter_4343" name="KmF6P" value="0.147"/>
          <Constant key="Parameter_4342" name="KmG6P" value="0.28"/>
          <Constant key="Parameter_4341" name="KmPEP" value="1.999"/>
          <Constant key="Parameter_4340" name="Vmax" value="2.32456"/>
          <Constant key="Parameter_4339" name="KmPGN" value="0.515958"/>
        </ListOfConstants>
        <KineticLaw function="Function_40">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_258">
              <SourceParameter reference="Metabolite_17"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_254">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_266">
              <SourceParameter reference="Parameter_4344"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_270">
              <SourceParameter reference="Parameter_4343"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_268">
              <SourceParameter reference="Parameter_4342"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_274">
              <SourceParameter reference="Parameter_4341"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_272">
              <SourceParameter reference="Metabolite_47"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_264">
              <SourceParameter reference="Parameter_4340"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_276">
              <SourceParameter reference="Metabolite_53"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_278">
              <SourceParameter reference="Parameter_4339"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_1" name="PFK" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_1">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-07T10:52:04Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_105" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_17" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_101" stoichiometry="1"/>
          <Product metabolite="Metabolite_19" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_47" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_117" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_119" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4338" name="KefrADP" value="0.0735264"/>
          <Constant key="Parameter_4337" name="KefrPEP" value="19.98"/>
          <Constant key="Parameter_4336" name="KeftADP" value="9.009"/>
          <Constant key="Parameter_4335" name="KeftPEP" value="0.26026"/>
          <Constant key="Parameter_4334" name="Keq" value="1998"/>
          <Constant key="Parameter_4333" name="KirADP" value="54.945"/>
          <Constant key="Parameter_4332" name="KirATP" value="2.4975e-005"/>
          <Constant key="Parameter_4331" name="KirF6P" value="1.84615"/>
          <Constant key="Parameter_4330" name="KirFDP" value="0.045954"/>
          <Constant key="Parameter_4329" name="KitADP" value="80.08"/>
          <Constant key="Parameter_4328" name="KitATP" value="0.014014"/>
          <Constant key="Parameter_4327" name="KitF6P" value="0.00856856"/>
          <Constant key="Parameter_4326" name="KitFDP" value="50.5505"/>
          <Constant key="Parameter_4325" name="KmrADP" value="0.690009"/>
          <Constant key="Parameter_4324" name="KmrATPMg" value="8.12187e-005"/>
          <Constant key="Parameter_4323" name="KmrF6P" value="2.05205e-005"/>
          <Constant key="Parameter_4322" name="KmrFDP" value="10.01"/>
          <Constant key="Parameter_4321" name="KmtADP" value="2.002"/>
          <Constant key="Parameter_4320" name="KmtATPMg" value="3.34334"/>
          <Constant key="Parameter_4319" name="KmtF6P" value="32.967"/>
          <Constant key="Parameter_4318" name="KmtFDP" value="9.99"/>
          <Constant key="Parameter_4317" name="L0" value="14.0851"/>
          <Constant key="Parameter_4316" name="Vmax" value="0.185253"/>
          <Constant key="Parameter_4315" name="Wr" value="0.0237041"/>
          <Constant key="Parameter_4314" name="Wt" value="0.146735"/>
          <Constant key="Parameter_4313" name="n" value="4"/>
        </ListOfConstants>
        <KineticLaw function="Function_41">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_298">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_267">
              <SourceParameter reference="Metabolite_17"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_280">
              <SourceParameter reference="Metabolite_19"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_330">
              <SourceParameter reference="Parameter_4338"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_332">
              <SourceParameter reference="Parameter_4337"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_324">
              <SourceParameter reference="Parameter_4336"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_328">
              <SourceParameter reference="Parameter_4335"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_282">
              <SourceParameter reference="Parameter_4334"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_296">
              <SourceParameter reference="Parameter_4333"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_300">
              <SourceParameter reference="Parameter_4332"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_284">
              <SourceParameter reference="Parameter_4331"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_290">
              <SourceParameter reference="Parameter_4330"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_318">
              <SourceParameter reference="Parameter_4329"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_320">
              <SourceParameter reference="Parameter_4328"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_314">
              <SourceParameter reference="Parameter_4327"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_308">
              <SourceParameter reference="Parameter_4326"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_292">
              <SourceParameter reference="Parameter_4325"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_286">
              <SourceParameter reference="Parameter_4324"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_294">
              <SourceParameter reference="Parameter_4323"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_288">
              <SourceParameter reference="Parameter_4322"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_310">
              <SourceParameter reference="Parameter_4321"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_316">
              <SourceParameter reference="Parameter_4320"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_312">
              <SourceParameter reference="Parameter_4319"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_306">
              <SourceParameter reference="Parameter_4318"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_304">
              <SourceParameter reference="Parameter_4317"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_246">
              <SourceParameter reference="Metabolite_117"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_271">
              <SourceParameter reference="Metabolite_119"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_326">
              <SourceParameter reference="Metabolite_47"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_279">
              <SourceParameter reference="Parameter_4316"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_302">
              <SourceParameter reference="Parameter_4315"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_322">
              <SourceParameter reference="Parameter_4314"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_275">
              <SourceParameter reference="Parameter_4313"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_2" name="FBA" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_2">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-07T10:49:12Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_19" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_13" stoichiometry="1"/>
          <Product metabolite="Metabolite_25" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_47" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4312" name="Keq" value="0.18981"/>
          <Constant key="Parameter_4311" name="KmDAP" value="0.13001"/>
          <Constant key="Parameter_4310" name="KmFDP" value="0.12012"/>
          <Constant key="Parameter_4309" name="KmGAP" value="0.13001"/>
          <Constant key="Parameter_4308" name="KmPEP" value="0.5"/>
          <Constant key="Parameter_4307" name="Vmax" value="21.6978"/>
        </ListOfConstants>
        <KineticLaw function="Function_42">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_325">
              <SourceParameter reference="Metabolite_13"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_329">
              <SourceParameter reference="Metabolite_19"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_321">
              <SourceParameter reference="Metabolite_25"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_317">
              <SourceParameter reference="Parameter_4312"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_309">
              <SourceParameter reference="Parameter_4311"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_313">
              <SourceParameter reference="Parameter_4310"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_305">
              <SourceParameter reference="Parameter_4309"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_297">
              <SourceParameter reference="Parameter_4308"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_301">
              <SourceParameter reference="Metabolite_47"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_333">
              <SourceParameter reference="Parameter_4307"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_3" name="TPI" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_3">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-07T10:49:41Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_13" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_25" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4306" name="Keq" value="0.270203"/>
          <Constant key="Parameter_4305" name="KmDAP" value="0.01"/>
          <Constant key="Parameter_4304" name="KmGAP" value="1.89301"/>
          <Constant key="Parameter_4303" name="Vmax" value="24.1843"/>
        </ListOfConstants>
        <KineticLaw function="Function_43">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_303">
              <SourceParameter reference="Metabolite_13"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_311">
              <SourceParameter reference="Metabolite_25"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_319">
              <SourceParameter reference="Parameter_4306"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_327">
              <SourceParameter reference="Parameter_4305"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_293">
              <SourceParameter reference="Parameter_4304"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_295">
              <SourceParameter reference="Parameter_4303"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_4" name="GDH" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_4">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T17:03:28Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_25" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_37" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_113" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_9" stoichiometry="1"/>
          <Product metabolite="Metabolite_39" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4302" name="Keq" value="20"/>
          <Constant key="Parameter_4301" name="KmBPG" value="0.2"/>
          <Constant key="Parameter_4300" name="KmGAP" value="2.47265"/>
          <Constant key="Parameter_4299" name="KmNAD" value="0.0110454"/>
          <Constant key="Parameter_4298" name="KmNADH" value="3.69797"/>
          <Constant key="Parameter_4297" name="KmP" value="0.017"/>
          <Constant key="Parameter_4296" name="Vmax" value="8.66573"/>
        </ListOfConstants>
        <KineticLaw function="Function_44">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_285">
              <SourceParameter reference="Metabolite_9"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_307">
              <SourceParameter reference="Metabolite_25"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_265">
              <SourceParameter reference="Parameter_4302"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_338">
              <SourceParameter reference="Parameter_4301"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_334">
              <SourceParameter reference="Parameter_4300"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_336">
              <SourceParameter reference="Parameter_4299"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_340">
              <SourceParameter reference="Parameter_4298"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_273">
              <SourceParameter reference="Parameter_4297"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_289">
              <SourceParameter reference="Metabolite_37"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_281">
              <SourceParameter reference="Metabolite_39"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_323">
              <SourceParameter reference="Metabolite_113"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_291">
              <SourceParameter reference="Parameter_4296"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_5" name="PGK" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_5">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-14T17:46:47Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_101" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_9" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_105" stoichiometry="1"/>
          <Product metabolite="Metabolite_51" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_117" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_119" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4295" name="Keq" value="99.9925"/>
          <Constant key="Parameter_4294" name="KmADPMg" value="0.085416"/>
          <Constant key="Parameter_4293" name="KmATPMg" value="3.47737"/>
          <Constant key="Parameter_4292" name="KmBPG" value="0.0113296"/>
          <Constant key="Parameter_4291" name="KmPGA3" value="2.45722"/>
          <Constant key="Parameter_4290" name="Vmax" value="16.1089"/>
        </ListOfConstants>
        <KineticLaw function="Function_45">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_277">
              <SourceParameter reference="Metabolite_9"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_315">
              <SourceParameter reference="Parameter_4295"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_342">
              <SourceParameter reference="Parameter_4294"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_346">
              <SourceParameter reference="Parameter_4293"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_344">
              <SourceParameter reference="Parameter_4292"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_348">
              <SourceParameter reference="Parameter_4291"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_337">
              <SourceParameter reference="Metabolite_117"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_262">
              <SourceParameter reference="Metabolite_119"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_287">
              <SourceParameter reference="Metabolite_51"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_341">
              <SourceParameter reference="Parameter_4290"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_6" name="GPM" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_6">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-24T12:56:58Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_51" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_49" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4289" name="Keq" value="0.565818"/>
          <Constant key="Parameter_4288" name="KmPGA2" value="1.9153"/>
          <Constant key="Parameter_4287" name="KmPGA3" value="0.115"/>
          <Constant key="Parameter_4286" name="Vmax" value="10.9934"/>
        </ListOfConstants>
        <KineticLaw function="Function_46">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_283">
              <SourceParameter reference="Parameter_4289"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_350">
              <SourceParameter reference="Parameter_4288"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_335">
              <SourceParameter reference="Parameter_4287"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_331">
              <SourceParameter reference="Metabolite_49"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_345">
              <SourceParameter reference="Metabolite_51"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_349">
              <SourceParameter reference="Parameter_4286"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_7" name="ENO" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_7">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T20:00:11Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_49" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_47" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4285" name="Keq" value="3"/>
          <Constant key="Parameter_4284" name="KmPEP" value="0.1"/>
          <Constant key="Parameter_4283" name="KmPGA2" value="0.1"/>
          <Constant key="Parameter_4282" name="Vmax" value="11.7189"/>
        </ListOfConstants>
        <KineticLaw function="Function_47">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_352">
              <SourceParameter reference="Parameter_4285"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_356">
              <SourceParameter reference="Parameter_4284"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_354">
              <SourceParameter reference="Parameter_4283"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_343">
              <SourceParameter reference="Metabolite_47"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_269">
              <SourceParameter reference="Metabolite_49"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_351">
              <SourceParameter reference="Parameter_4282"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_8" name="PYK" reversible="false" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_8">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-08T11:11:05Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_101" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_47" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_105" stoichiometry="1"/>
          <Product metabolite="Metabolite_55" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_19" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_23" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_27" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_117" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_119" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_61" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_63" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_65" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_69" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_71" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4281" name="KefrFDP" value="0.449149"/>
          <Constant key="Parameter_4280" name="KefrG6P" value="0.158746"/>
          <Constant key="Parameter_4279" name="KefrGL6P" value="0.150482"/>
          <Constant key="Parameter_4278" name="KefrR5P" value="9.33254"/>
          <Constant key="Parameter_4277" name="KefrRU5P" value="1.53591"/>
          <Constant key="Parameter_4276" name="KefrS7P" value="0.0785955"/>
          <Constant key="Parameter_4275" name="KefrX5P" value="0.677374"/>
          <Constant key="Parameter_4274" name="KeftATP" value="3.69117"/>
          <Constant key="Parameter_4273" name="KeftSUCCOA" value="8.26406"/>
          <Constant key="Parameter_4272" name="KirADP" value="0.517585"/>
          <Constant key="Parameter_4271" name="KirATP" value="96.0333"/>
          <Constant key="Parameter_4270" name="KirPEP" value="0.181056"/>
          <Constant key="Parameter_4269" name="KirPYR" value="15.1403"/>
          <Constant key="Parameter_4268" name="KirPyrATP" value="230.781"/>
          <Constant key="Parameter_4267" name="KitADP" value="0.224911"/>
          <Constant key="Parameter_4266" name="KitATP" value="0.039564"/>
          <Constant key="Parameter_4265" name="KitPEP" value="0.465672"/>
          <Constant key="Parameter_4264" name="KitPYR" value="0.2499"/>
          <Constant key="Parameter_4263" name="KitPyrATP" value="11.3691"/>
          <Constant key="Parameter_4262" name="KmrADPMg" value="0.326144"/>
          <Constant key="Parameter_4261" name="KmrPEP" value="5.56368e-007"/>
          <Constant key="Parameter_4260" name="KmtADPMg" value="0.054678"/>
          <Constant key="Parameter_4259" name="KmtPEP" value="0.11475"/>
          <Constant key="Parameter_4258" name="L0" value="50.4818"/>
          <Constant key="Parameter_4257" name="Vmax" value="0.74716"/>
          <Constant key="Parameter_4256" name="n" value="4"/>
        </ListOfConstants>
        <KineticLaw function="Function_48">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_370">
              <SourceParameter reference="Metabolite_101"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_402">
              <SourceParameter reference="Metabolite_19"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_406">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_410">
              <SourceParameter reference="Metabolite_27"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_404">
              <SourceParameter reference="Parameter_4281"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_408">
              <SourceParameter reference="Parameter_4280"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_412">
              <SourceParameter reference="Parameter_4279"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_416">
              <SourceParameter reference="Parameter_4278"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_420">
              <SourceParameter reference="Parameter_4277"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_424">
              <SourceParameter reference="Parameter_4276"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_428">
              <SourceParameter reference="Parameter_4275"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_400">
              <SourceParameter reference="Parameter_4274"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_398">
              <SourceParameter reference="Parameter_4273"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_372">
              <SourceParameter reference="Parameter_4272"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_368">
              <SourceParameter reference="Parameter_4271"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_360">
              <SourceParameter reference="Parameter_4270"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_376">
              <SourceParameter reference="Parameter_4269"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_378">
              <SourceParameter reference="Parameter_4268"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_390">
              <SourceParameter reference="Parameter_4267"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_388">
              <SourceParameter reference="Parameter_4266"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_384">
              <SourceParameter reference="Parameter_4265"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_392">
              <SourceParameter reference="Parameter_4264"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_394">
              <SourceParameter reference="Parameter_4263"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_362">
              <SourceParameter reference="Parameter_4262"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_364">
              <SourceParameter reference="Parameter_4261"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_386">
              <SourceParameter reference="Parameter_4260"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_382">
              <SourceParameter reference="Parameter_4259"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_380">
              <SourceParameter reference="Parameter_4258"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_358">
              <SourceParameter reference="Metabolite_117"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_366">
              <SourceParameter reference="Metabolite_119"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_299">
              <SourceParameter reference="Metabolite_47"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_374">
              <SourceParameter reference="Metabolite_55"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_414">
              <SourceParameter reference="Metabolite_61"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_418">
              <SourceParameter reference="Metabolite_63"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_422">
              <SourceParameter reference="Metabolite_65"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_396">
              <SourceParameter reference="Metabolite_69"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_357">
              <SourceParameter reference="Parameter_4257"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_426">
              <SourceParameter reference="Metabolite_71"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_353">
              <SourceParameter reference="Parameter_4256"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_9" name="ZWF" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_9">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:06:06Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_23" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_41" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_27" stoichiometry="1"/>
          <Product metabolite="Metabolite_43" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4255" name="KdG6P" value="0.192"/>
          <Constant key="Parameter_4254" name="KdGL6P" value="0.02"/>
          <Constant key="Parameter_4253" name="Keq" value="6e+010"/>
          <Constant key="Parameter_4252" name="KmG6P" value="0.118525"/>
          <Constant key="Parameter_4251" name="KmGL6P" value="0.328629"/>
          <Constant key="Parameter_4250" name="KmNADP" value="0.0274"/>
          <Constant key="Parameter_4249" name="KmNADPH" value="0.0168"/>
          <Constant key="Parameter_4248" name="Vmax" value="0.2658"/>
        </ListOfConstants>
        <KineticLaw function="Function_49">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_425">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_417">
              <SourceParameter reference="Metabolite_27"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_405">
              <SourceParameter reference="Parameter_4255"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_389">
              <SourceParameter reference="Parameter_4254"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_409">
              <SourceParameter reference="Parameter_4253"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_397">
              <SourceParameter reference="Parameter_4252"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_393">
              <SourceParameter reference="Parameter_4251"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_401">
              <SourceParameter reference="Parameter_4250"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_385">
              <SourceParameter reference="Parameter_4249"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_421">
              <SourceParameter reference="Metabolite_41"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_413">
              <SourceParameter reference="Metabolite_43"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_429">
              <SourceParameter reference="Parameter_4248"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_10" name="PGL" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_10">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-05T14:31:20Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_27" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_53" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_23" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4247" name="Keq" value="42.7572"/>
          <Constant key="Parameter_4246" name="KiG6P" value="2.0001"/>
          <Constant key="Parameter_4245" name="KmGL6P" value="0.022977"/>
          <Constant key="Parameter_4244" name="KmPGN" value="9.99"/>
          <Constant key="Parameter_4243" name="Vmax" value="11.5967"/>
        </ListOfConstants>
        <KineticLaw function="Function_50">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_381">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_391">
              <SourceParameter reference="Metabolite_27"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_407">
              <SourceParameter reference="Parameter_4247"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_377">
              <SourceParameter reference="Parameter_4246"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_415">
              <SourceParameter reference="Parameter_4245"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_423">
              <SourceParameter reference="Parameter_4244"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_399">
              <SourceParameter reference="Metabolite_53"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_383">
              <SourceParameter reference="Parameter_4243"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_11" name="GND" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_11">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:44:26Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_41" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_53" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_43" stoichiometry="1"/>
          <Product metabolite="Metabolite_63" stoichiometry="1"/>
          <Product metabolite="Metabolite_111" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_105" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_19" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4242" name="KdHCO3" value="58.8951"/>
          <Constant key="Parameter_4241" name="KdHCO3NADPH" value="9.7257"/>
          <Constant key="Parameter_4240" name="KdNADP" value="0.116989"/>
          <Constant key="Parameter_4239" name="KdNADPH" value="0.00340034"/>
          <Constant key="Parameter_4238" name="KdRu5P" value="0.0440044"/>
          <Constant key="Parameter_4237" name="KefATP" value="0.0650065"/>
          <Constant key="Parameter_4236" name="KefFbP" value="0.0129987"/>
          <Constant key="Parameter_4235" name="KefNADPATP" value="0.139986"/>
          <Constant key="Parameter_4234" name="KefNADPFbP" value="0.00519948"/>
          <Constant key="Parameter_4233" name="Keq" value="49.962"/>
          <Constant key="Parameter_4232" name="KmHCO3" value="6.41899"/>
          <Constant key="Parameter_4231" name="KmNADP" value="0.049"/>
          <Constant key="Parameter_4230" name="KmNADPH" value="68.3828"/>
          <Constant key="Parameter_4229" name="KmPGN" value="0.093"/>
          <Constant key="Parameter_4228" name="KmRU5P" value="45.1977"/>
          <Constant key="Parameter_4227" name="Vmax" value="4.08105"/>
        </ListOfConstants>
        <KineticLaw function="Function_51">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_437">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_355">
              <SourceParameter reference="Metabolite_19"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_369">
              <SourceParameter reference="Metabolite_111"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_443">
              <SourceParameter reference="Parameter_4242"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_449">
              <SourceParameter reference="Parameter_4241"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_361">
              <SourceParameter reference="Parameter_4240"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_445">
              <SourceParameter reference="Parameter_4239"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_447">
              <SourceParameter reference="Parameter_4238"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_439">
              <SourceParameter reference="Parameter_4237"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_431">
              <SourceParameter reference="Parameter_4236"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_441">
              <SourceParameter reference="Parameter_4235"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_433">
              <SourceParameter reference="Parameter_4234"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_365">
              <SourceParameter reference="Parameter_4233"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_455">
              <SourceParameter reference="Parameter_4232"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_435">
              <SourceParameter reference="Parameter_4231"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_451">
              <SourceParameter reference="Parameter_4230"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_339">
              <SourceParameter reference="Parameter_4229"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_453">
              <SourceParameter reference="Parameter_4228"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_427">
              <SourceParameter reference="Metabolite_41"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_395">
              <SourceParameter reference="Metabolite_43"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_411">
              <SourceParameter reference="Metabolite_53"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_373">
              <SourceParameter reference="Metabolite_63"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_375">
              <SourceParameter reference="Parameter_4227"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_12" name="RPE" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_12">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-21T14:42:32Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_63" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_71" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4226" name="Keq" value="1.5015"/>
          <Constant key="Parameter_4225" name="KmRU5P" value="0.872522"/>
          <Constant key="Parameter_4224" name="KmX5P" value="0.893607"/>
          <Constant key="Parameter_4223" name="Vmax" value="6.00103"/>
        </ListOfConstants>
        <KineticLaw function="Function_52">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_444">
              <SourceParameter reference="Parameter_4226"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_440">
              <SourceParameter reference="Parameter_4225"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_436">
              <SourceParameter reference="Parameter_4224"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_452">
              <SourceParameter reference="Metabolite_63"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_456">
              <SourceParameter reference="Parameter_4223"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_448">
              <SourceParameter reference="Metabolite_71"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_13" name="RPI" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_13">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-22T20:48:16Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_63" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_61" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_15" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4222" name="Keq" value="0.330093"/>
          <Constant key="Parameter_4221" name="KmE4P" value="0.67067"/>
          <Constant key="Parameter_4220" name="KmR5P" value="3.09715"/>
          <Constant key="Parameter_4219" name="KmRU5P" value="4.40263"/>
          <Constant key="Parameter_4218" name="Vmax" value="8"/>
        </ListOfConstants>
        <KineticLaw function="Function_53">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_371">
              <SourceParameter reference="Metabolite_15"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_432">
              <SourceParameter reference="Parameter_4222"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_403">
              <SourceParameter reference="Parameter_4221"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_363">
              <SourceParameter reference="Parameter_4220"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_347">
              <SourceParameter reference="Parameter_4219"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_450">
              <SourceParameter reference="Metabolite_61"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_442">
              <SourceParameter reference="Metabolite_63"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_434">
              <SourceParameter reference="Parameter_4218"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_14" name="X5P_GAP_TKT" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_14">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T11:51:30Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_97" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_71" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_25" stoichiometry="1"/>
          <Product metabolite="Metabolite_99" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4217" name="Keq" value="1.001"/>
          <Constant key="Parameter_4216" name="kcat" value="40"/>
        </ListOfConstants>
        <KineticLaw function="Function_54">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_446">
              <SourceParameter reference="Metabolite_25"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_458">
              <SourceParameter reference="Parameter_4217"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_430">
              <SourceParameter reference="Metabolite_71"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_419">
              <SourceParameter reference="Parameter_4216"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_367">
              <SourceParameter reference="Metabolite_97"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_379">
              <SourceParameter reference="Metabolite_99"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_15" name="F6P_E4P_TKT" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_15">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T20:00:11Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_15" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_99" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_17" stoichiometry="1"/>
          <Product metabolite="Metabolite_97" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4215" name="Keq" value="0.5005"/>
          <Constant key="Parameter_4214" name="kcat" value="40.0002"/>
        </ListOfConstants>
        <KineticLaw function="Function_55">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_438">
              <SourceParameter reference="Metabolite_15"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_460">
              <SourceParameter reference="Metabolite_17"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_464">
              <SourceParameter reference="Parameter_4215"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_459">
              <SourceParameter reference="Parameter_4214"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_462">
              <SourceParameter reference="Metabolite_97"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_359">
              <SourceParameter reference="Metabolite_99"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_16" name="S7P_R5P_TKT" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_16">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-09T18:41:37Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_61" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_99" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_65" stoichiometry="1"/>
          <Product metabolite="Metabolite_97" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4213" name="Keq" value="0.33033"/>
          <Constant key="Parameter_4212" name="kcat" value="199.979"/>
        </ListOfConstants>
        <KineticLaw function="Function_56">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_470">
              <SourceParameter reference="Parameter_4213"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_461">
              <SourceParameter reference="Metabolite_61"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_466">
              <SourceParameter reference="Metabolite_65"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_465">
              <SourceParameter reference="Parameter_4212"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_468">
              <SourceParameter reference="Metabolite_97"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_454">
              <SourceParameter reference="Metabolite_99"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_17" name="F6P_GAP_TAL" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_17">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-27T11:53:43Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_25" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_95" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_17" stoichiometry="1"/>
          <Product metabolite="Metabolite_93" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4211" name="Keq" value="0.11011"/>
          <Constant key="Parameter_4210" name="kcat" value="119.992"/>
        </ListOfConstants>
        <KineticLaw function="Function_57">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_472">
              <SourceParameter reference="Metabolite_17"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_467">
              <SourceParameter reference="Metabolite_25"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_476">
              <SourceParameter reference="Parameter_4211"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_471">
              <SourceParameter reference="Parameter_4210"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_474">
              <SourceParameter reference="Metabolite_93"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_387">
              <SourceParameter reference="Metabolite_95"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_18" name="S7P_E4P_TAL" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_18">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-09T18:41:34Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_65" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_93" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_15" stoichiometry="1"/>
          <Product metabolite="Metabolite_95" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4209" name="Keq" value="26.6266"/>
          <Constant key="Parameter_4208" name="kcat" value="99.9991"/>
        </ListOfConstants>
        <KineticLaw function="Function_58">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_478">
              <SourceParameter reference="Metabolite_15"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_482">
              <SourceParameter reference="Parameter_4209"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_473">
              <SourceParameter reference="Metabolite_65"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_477">
              <SourceParameter reference="Parameter_4208"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_457">
              <SourceParameter reference="Metabolite_93"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_480">
              <SourceParameter reference="Metabolite_95"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_19" name="FBP" reversible="false" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_19">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T17:07:07Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_19" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_17" stoichiometry="1"/>
          <Product metabolite="Metabolite_113" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_103" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_115" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_121" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4207" name="KirAMP" value="0.00122122"/>
          <Constant key="Parameter_4206" name="KirAMPFDP" value="0.256256"/>
          <Constant key="Parameter_4205" name="KirF6P" value="1.12112"/>
          <Constant key="Parameter_4204" name="KirF6PMg" value="0.384615"/>
          <Constant key="Parameter_4203" name="KirFDP" value="1.35327"/>
          <Constant key="Parameter_4202" name="KirFDPMg" value="0.75924"/>
          <Constant key="Parameter_4201" name="KirFDPMgMg" value="0.356356"/>
          <Constant key="Parameter_4200" name="KirP" value="3.16316"/>
          <Constant key="Parameter_4199" name="KirPF6P" value="6.60538"/>
          <Constant key="Parameter_4198" name="KirPF6PMg" value="48.4484"/>
          <Constant key="Parameter_4197" name="KirPMg" value="0.856"/>
          <Constant key="Parameter_4196" name="KitAMP" value="0.000255"/>
          <Constant key="Parameter_4195" name="KitAMPFDP" value="690"/>
          <Constant key="Parameter_4194" name="KitF6P" value="0.304"/>
          <Constant key="Parameter_4193" name="KitF6PMg" value="315"/>
          <Constant key="Parameter_4192" name="KitFDP" value="0.043101"/>
          <Constant key="Parameter_4191" name="KitFDPMg" value="0.00642"/>
          <Constant key="Parameter_4190" name="KitFDPMgMg" value="100"/>
          <Constant key="Parameter_4189" name="KitP" value="0.642"/>
          <Constant key="Parameter_4188" name="KitPF6P" value="0.00689"/>
          <Constant key="Parameter_4187" name="KitPF6PMg" value="16.5"/>
          <Constant key="Parameter_4186" name="KitPMg" value="539"/>
          <Constant key="Parameter_4185" name="KmrFDP" value="0.0636141"/>
          <Constant key="Parameter_4184" name="KmrMg" value="0.039039"/>
          <Constant key="Parameter_4183" name="KmtFDP" value="1e-005"/>
          <Constant key="Parameter_4182" name="KmtMg" value="55.055"/>
          <Constant key="Parameter_4181" name="L0" value="0.000815"/>
          <Constant key="Parameter_4180" name="Vmax" value="0.215583"/>
          <Constant key="Parameter_4179" name="n" value="4"/>
          <Constant key="Parameter_4178" name="KdFDPMg" value="5.81"/>
        </ListOfConstants>
        <KineticLaw function="Function_59">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_514">
              <SourceParameter reference="Metabolite_103"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_500">
              <SourceParameter reference="Metabolite_17"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_510">
              <SourceParameter reference="Metabolite_19"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_512">
              <SourceParameter reference="ModelValue_3"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_516">
              <SourceParameter reference="Parameter_4207"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_520">
              <SourceParameter reference="Parameter_4206"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_502">
              <SourceParameter reference="Parameter_4205"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_504">
              <SourceParameter reference="Parameter_4204"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_488">
              <SourceParameter reference="Parameter_4203"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_484">
              <SourceParameter reference="Parameter_4202"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_518">
              <SourceParameter reference="Parameter_4201"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_496">
              <SourceParameter reference="Parameter_4200"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_506">
              <SourceParameter reference="Parameter_4199"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_508">
              <SourceParameter reference="Parameter_4198"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_498">
              <SourceParameter reference="Parameter_4197"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_542">
              <SourceParameter reference="Parameter_4196"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_548">
              <SourceParameter reference="Parameter_4195"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_534">
              <SourceParameter reference="Parameter_4194"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_536">
              <SourceParameter reference="Parameter_4193"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_526">
              <SourceParameter reference="Parameter_4192"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_544">
              <SourceParameter reference="Parameter_4191"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_546">
              <SourceParameter reference="Parameter_4190"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_530">
              <SourceParameter reference="Parameter_4189"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_538">
              <SourceParameter reference="Parameter_4188"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_540">
              <SourceParameter reference="Parameter_4187"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_532">
              <SourceParameter reference="Parameter_4186"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_486">
              <SourceParameter reference="Parameter_4185"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_492">
              <SourceParameter reference="Parameter_4184"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_524">
              <SourceParameter reference="Parameter_4183"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_528">
              <SourceParameter reference="Parameter_4182"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_522">
              <SourceParameter reference="Parameter_4181"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_490">
              <SourceParameter reference="Metabolite_115"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_463">
              <SourceParameter reference="Metabolite_121"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_494">
              <SourceParameter reference="Metabolite_113"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_483">
              <SourceParameter reference="Parameter_4180"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_479">
              <SourceParameter reference="Parameter_4179"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_20" name="PPC" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_20">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:45:27Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_47" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_111" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_45" stoichiometry="1"/>
          <Product metabolite="Metabolite_113" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_1" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_11" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_19" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_21" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_35" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_67" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_123" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_125" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4177" name="KdrOAA" value="4.35404"/>
          <Constant key="Parameter_4176" name="KdrPEP" value="655.591"/>
          <Constant key="Parameter_4175" name="KdtOAA" value="17.9127"/>
          <Constant key="Parameter_4174" name="KdtPEP" value="0.0121991"/>
          <Constant key="Parameter_4173" name="KefrACCOA" value="0.140029"/>
          <Constant key="Parameter_4172" name="KefrASP" value="0.389183"/>
          <Constant key="Parameter_4171" name="KefrCIT" value="34.4277"/>
          <Constant key="Parameter_4170" name="KefrCYS" value="0.000448843"/>
          <Constant key="Parameter_4169" name="KefrFDP" value="9.99126"/>
          <Constant key="Parameter_4168" name="KefrFDPACCOA" value="0.0156251"/>
          <Constant key="Parameter_4167" name="KefrFUM" value="2.7475"/>
          <Constant key="Parameter_4166" name="KefrMAL" value="0.23002"/>
          <Constant key="Parameter_4165" name="KefrSUC" value="22.9834"/>
          <Constant key="Parameter_4164" name="KeftACCOA" value="1.27598"/>
          <Constant key="Parameter_4163" name="KeftASP" value="27.4729"/>
          <Constant key="Parameter_4162" name="KeftCIT" value="0.521945"/>
          <Constant key="Parameter_4161" name="KeftCYS" value="0.977374"/>
          <Constant key="Parameter_4160" name="KeftFDP" value="13.2064"/>
          <Constant key="Parameter_4159" name="KeftFDPACCOA" value="47.7563"/>
          <Constant key="Parameter_4158" name="KeftFUM" value="9.76119"/>
          <Constant key="Parameter_4157" name="KeftMAL" value="0.737283"/>
          <Constant key="Parameter_4156" name="KeftSUC" value="107.18"/>
          <Constant key="Parameter_4155" name="Keq" value="149.705"/>
          <Constant key="Parameter_4154" name="KmrHCO3" value="0.00219811"/>
          <Constant key="Parameter_4153" name="KmrOAA" value="13.0469"/>
          <Constant key="Parameter_4152" name="KmrP" value="0.663356"/>
          <Constant key="Parameter_4151" name="KmrPEP" value="3.20089"/>
          <Constant key="Parameter_4150" name="KmtHCO3" value="0.00220031"/>
          <Constant key="Parameter_4149" name="KmtOAA" value="6.80995"/>
          <Constant key="Parameter_4148" name="KmtP" value="0.285131"/>
          <Constant key="Parameter_4147" name="KmtPEP" value="5.12497"/>
          <Constant key="Parameter_4146" name="L0" value="6.37209e-006"/>
          <Constant key="Parameter_4145" name="Vmax" value="21.439"/>
          <Constant key="Parameter_4144" name="n" value="4.00319"/>
        </ListOfConstants>
        <KineticLaw function="Function_60">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_493">
              <SourceParameter reference="Metabolite_1"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_564">
              <SourceParameter reference="Metabolite_123"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_572">
              <SourceParameter reference="Metabolite_11"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_568">
              <SourceParameter reference="Metabolite_125"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_485">
              <SourceParameter reference="Metabolite_19"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_576">
              <SourceParameter reference="Metabolite_21"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_537">
              <SourceParameter reference="Metabolite_111"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_505">
              <SourceParameter reference="Parameter_4177"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_521">
              <SourceParameter reference="Parameter_4176"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_560">
              <SourceParameter reference="Parameter_4175"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_554">
              <SourceParameter reference="Parameter_4174"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_588">
              <SourceParameter reference="Parameter_4173"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_594">
              <SourceParameter reference="Parameter_4172"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_598">
              <SourceParameter reference="Parameter_4171"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_596">
              <SourceParameter reference="Parameter_4170"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_590">
              <SourceParameter reference="Parameter_4169"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_592">
              <SourceParameter reference="Parameter_4168"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_600">
              <SourceParameter reference="Parameter_4167"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_602">
              <SourceParameter reference="Parameter_4166"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_604">
              <SourceParameter reference="Parameter_4165"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_489">
              <SourceParameter reference="Parameter_4164"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_566">
              <SourceParameter reference="Parameter_4163"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_574">
              <SourceParameter reference="Parameter_4162"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_570">
              <SourceParameter reference="Parameter_4161"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_469">
              <SourceParameter reference="Parameter_4160"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_550">
              <SourceParameter reference="Parameter_4159"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_578">
              <SourceParameter reference="Parameter_4158"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_582">
              <SourceParameter reference="Parameter_4157"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_586">
              <SourceParameter reference="Parameter_4156"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_525">
              <SourceParameter reference="Parameter_4155"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_517">
              <SourceParameter reference="Parameter_4154"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_509">
              <SourceParameter reference="Parameter_4153"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_501">
              <SourceParameter reference="Parameter_4152"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_513">
              <SourceParameter reference="Parameter_4151"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_556">
              <SourceParameter reference="Parameter_4150"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_558">
              <SourceParameter reference="Parameter_4149"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_562">
              <SourceParameter reference="Parameter_4148"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_552">
              <SourceParameter reference="Parameter_4147"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_497">
              <SourceParameter reference="Parameter_4146"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_580">
              <SourceParameter reference="Metabolite_35"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_533">
              <SourceParameter reference="Metabolite_45"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_529">
              <SourceParameter reference="Metabolite_113"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_541">
              <SourceParameter reference="Metabolite_47"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_584">
              <SourceParameter reference="Metabolite_67"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_549">
              <SourceParameter reference="Parameter_4145"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_545">
              <SourceParameter reference="Parameter_4144"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_21" name="PCK" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_21">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:46:44Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_105" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_45" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_101" stoichiometry="1"/>
          <Product metabolite="Metabolite_47" stoichiometry="1"/>
          <Product metabolite="Metabolite_111" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_117" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_119" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4143" name="Keq" value="1.88166"/>
          <Constant key="Parameter_4142" name="KmADP" value="0.04997"/>
          <Constant key="Parameter_4141" name="KmATP" value="0.0600418"/>
          <Constant key="Parameter_4140" name="KmHCO3" value="2.6319"/>
          <Constant key="Parameter_4139" name="KmOAA" value="0.66966"/>
          <Constant key="Parameter_4138" name="KmPEP" value="0.0699833"/>
          <Constant key="Parameter_4137" name="Vmax" value="8.08777"/>
        </ListOfConstants>
        <KineticLaw function="Function_61">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_565">
              <SourceParameter reference="Metabolite_101"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_593">
              <SourceParameter reference="Metabolite_111"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_581">
              <SourceParameter reference="Parameter_4143"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_561">
              <SourceParameter reference="Parameter_4142"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_577">
              <SourceParameter reference="Parameter_4141"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_569">
              <SourceParameter reference="Parameter_4140"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_573">
              <SourceParameter reference="Parameter_4139"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_557">
              <SourceParameter reference="Parameter_4138"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_589">
              <SourceParameter reference="Metabolite_117"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_601">
              <SourceParameter reference="Metabolite_119"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_597">
              <SourceParameter reference="Metabolite_45"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_585">
              <SourceParameter reference="Metabolite_47"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_605">
              <SourceParameter reference="Parameter_4137"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_22" name="PPS" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_22">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:57:27Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_105" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_55" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_103" stoichiometry="1"/>
          <Product metabolite="Metabolite_47" stoichiometry="1"/>
          <Product metabolite="Metabolite_113" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_101" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_7" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_45" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_115" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_117" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_119" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4136" name="KdAMP" value="1480"/>
          <Constant key="Parameter_4135" name="KdATPMgPPS" value="0.0549"/>
          <Constant key="Parameter_4134" name="KdMg" value="36.9"/>
          <Constant key="Parameter_4133" name="KdP" value="346"/>
          <Constant key="Parameter_4132" name="KdPEP" value="95.7"/>
          <Constant key="Parameter_4131" name="KdPYR" value="2740"/>
          <Constant key="Parameter_4130" name="KefADP" value="0.0283"/>
          <Constant key="Parameter_4129" name="KefAKG" value="0.274"/>
          <Constant key="Parameter_4128" name="KefATP" value="0.000628"/>
          <Constant key="Parameter_4127" name="KefOAA" value="0.796"/>
          <Constant key="Parameter_4126" name="Keq" value="200000"/>
          <Constant key="Parameter_4125" name="KmAMP" value="0.000384"/>
          <Constant key="Parameter_4124" name="KmATPMg" value="0.0549"/>
          <Constant key="Parameter_4123" name="KmP" value="84.4"/>
          <Constant key="Parameter_4122" name="KmPEP" value="20.7"/>
          <Constant key="Parameter_4121" name="KmPYR" value="0.229"/>
          <Constant key="Parameter_4120" name="Vmax" value="0.0163772"/>
          <Constant key="Parameter_4119" name="W" value="10"/>
          <Constant key="Parameter_4118" name="alpha" value="38900"/>
          <Constant key="Parameter_4117" name="KdADPMg" value="1.27771"/>
          <Constant key="Parameter_4116" name="KdATPMg" value="0.0847634"/>
        </ListOfConstants>
        <KineticLaw function="Function_62">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_625">
              <SourceParameter reference="Metabolite_101"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_607">
              <SourceParameter reference="Metabolite_7"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_579">
              <SourceParameter reference="Metabolite_103"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_631">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_619">
              <SourceParameter reference="ModelValue_1"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_515">
              <SourceParameter reference="Parameter_4136"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_627">
              <SourceParameter reference="ModelValue_2"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_539">
              <SourceParameter reference="Parameter_4135"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_523">
              <SourceParameter reference="Parameter_4134"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_507">
              <SourceParameter reference="Parameter_4133"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_633">
              <SourceParameter reference="Parameter_4132"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_635">
              <SourceParameter reference="Parameter_4131"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_623">
              <SourceParameter reference="Parameter_4130"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_609">
              <SourceParameter reference="Parameter_4129"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_629">
              <SourceParameter reference="Parameter_4128"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_615">
              <SourceParameter reference="Parameter_4127"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_551">
              <SourceParameter reference="Parameter_4126"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_617">
              <SourceParameter reference="Parameter_4125"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_475">
              <SourceParameter reference="Parameter_4124"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_531">
              <SourceParameter reference="Parameter_4123"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_611">
              <SourceParameter reference="Parameter_4122"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_491">
              <SourceParameter reference="Parameter_4121"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_603">
              <SourceParameter reference="Metabolite_115"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_621">
              <SourceParameter reference="Metabolite_117"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_563">
              <SourceParameter reference="Metabolite_119"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_613">
              <SourceParameter reference="Metabolite_45"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_595">
              <SourceParameter reference="Metabolite_113"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_587">
              <SourceParameter reference="Metabolite_47"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_571">
              <SourceParameter reference="Metabolite_55"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_555">
              <SourceParameter reference="Parameter_4120"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_547">
              <SourceParameter reference="Parameter_4119"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_499">
              <SourceParameter reference="Parameter_4118"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_23" name="MAE" reversible="false" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_23">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:56:51Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_35" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_37" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_39" stoichiometry="1"/>
          <Product metabolite="Metabolite_55" stoichiometry="1"/>
          <Product metabolite="Metabolite_111" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_105" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_1" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_109" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_115" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_123" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_127" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4115" name="KefrACCOA" value="1.83336"/>
          <Constant key="Parameter_4114" name="KefrASP" value="0.362142"/>
          <Constant key="Parameter_4113" name="KefrATP" value="88.9752"/>
          <Constant key="Parameter_4112" name="KefrCOA" value="0.268"/>
          <Constant key="Parameter_4111" name="KeftACCOA" value="0.197"/>
          <Constant key="Parameter_4110" name="KeftASP" value="0.583"/>
          <Constant key="Parameter_4109" name="KeftATP" value="0.26"/>
          <Constant key="Parameter_4108" name="KeftCOA" value="0.268"/>
          <Constant key="Parameter_4107" name="KirNAD" value="0.636457"/>
          <Constant key="Parameter_4106" name="KitNAD" value="0.990398"/>
          <Constant key="Parameter_4105" name="KmrMAL" value="0.212913"/>
          <Constant key="Parameter_4104" name="KmrMg" value="0.191871"/>
          <Constant key="Parameter_4103" name="KmrMn" value="0.272568"/>
          <Constant key="Parameter_4102" name="KmrNAD" value="1.36636"/>
          <Constant key="Parameter_4101" name="KmtMAL" value="0.093"/>
          <Constant key="Parameter_4100" name="KmtMg" value="2.37681"/>
          <Constant key="Parameter_4099" name="KmtMn" value="0.410198"/>
          <Constant key="Parameter_4098" name="KmtNAD" value="0.108"/>
          <Constant key="Parameter_4097" name="L0" value="19.9"/>
          <Constant key="Parameter_4096" name="Vmax" value="6.64269"/>
          <Constant key="Parameter_4095" name="n" value="4"/>
        </ListOfConstants>
        <KineticLaw function="Function_63">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_641">
              <SourceParameter reference="Metabolite_1"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_481">
              <SourceParameter reference="Metabolite_123"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_637">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_645">
              <SourceParameter reference="Metabolite_109"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_659">
              <SourceParameter reference="Parameter_4115"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_655">
              <SourceParameter reference="Parameter_4114"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_657">
              <SourceParameter reference="Parameter_4113"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_661">
              <SourceParameter reference="Parameter_4112"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_643">
              <SourceParameter reference="Parameter_4111"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_599">
              <SourceParameter reference="Parameter_4110"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_639">
              <SourceParameter reference="Parameter_4109"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_647">
              <SourceParameter reference="Parameter_4108"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_616">
              <SourceParameter reference="Parameter_4107"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_651">
              <SourceParameter reference="Parameter_4106"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_620">
              <SourceParameter reference="Parameter_4105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_543">
              <SourceParameter reference="Parameter_4104"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_511">
              <SourceParameter reference="Parameter_4103"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_612">
              <SourceParameter reference="Parameter_4102"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_653">
              <SourceParameter reference="Parameter_4101"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_583">
              <SourceParameter reference="Parameter_4100"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_567">
              <SourceParameter reference="Parameter_4099"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_649">
              <SourceParameter reference="Parameter_4098"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_495">
              <SourceParameter reference="Parameter_4097"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_628">
              <SourceParameter reference="Metabolite_35"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_608">
              <SourceParameter reference="Metabolite_115"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_527">
              <SourceParameter reference="Metabolite_127"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_624">
              <SourceParameter reference="Metabolite_37"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_636">
              <SourceParameter reference="Parameter_4096"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_632">
              <SourceParameter reference="Parameter_4095"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_24" name="PDH" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_24">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:47:07Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_109" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_37" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_55" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_1" stoichiometry="1"/>
          <Product metabolite="Metabolite_39" stoichiometry="1"/>
          <Product metabolite="Metabolite_111" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4094" name="Keq" value="3138.16"/>
          <Constant key="Parameter_4093" name="KmACCOA" value="10.174"/>
          <Constant key="Parameter_4092" name="KmCOA" value="0.00500461"/>
          <Constant key="Parameter_4091" name="KmHCO3" value="0.00545112"/>
          <Constant key="Parameter_4090" name="KmNAD" value="0.00999"/>
          <Constant key="Parameter_4089" name="KmNADH" value="6.63512"/>
          <Constant key="Parameter_4088" name="KmPYR" value="2"/>
          <Constant key="Parameter_4087" name="Vmax" value="961.706"/>
        </ListOfConstants>
        <KineticLaw function="Function_64">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_646">
              <SourceParameter reference="Metabolite_1"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_658">
              <SourceParameter reference="Metabolite_109"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_638">
              <SourceParameter reference="Metabolite_111"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_575">
              <SourceParameter reference="Parameter_4094"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_610">
              <SourceParameter reference="Parameter_4093"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_553">
              <SourceParameter reference="Parameter_4092"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_626">
              <SourceParameter reference="Parameter_4091"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_503">
              <SourceParameter reference="Parameter_4090"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_618">
              <SourceParameter reference="Parameter_4089"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_535">
              <SourceParameter reference="Parameter_4088"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_654">
              <SourceParameter reference="Metabolite_37"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_642">
              <SourceParameter reference="Metabolite_39"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_650">
              <SourceParameter reference="Metabolite_55"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_662">
              <SourceParameter reference="Parameter_4087"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_25" name="GLT" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_25">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:45:20Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_1" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_45" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_11" stoichiometry="1"/>
          <Product metabolite="Metabolite_109" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_105" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_7" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_39" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4086" name="KdACCOA0" value="0.699309"/>
          <Constant key="Parameter_4085" name="KdcsCIT" value="7.38128"/>
          <Constant key="Parameter_4084" name="KdcsCOA" value="0.00174282"/>
          <Constant key="Parameter_4083" name="KdcsOAA" value="0.154601"/>
          <Constant key="Parameter_4082" name="Keq" value="83067.1"/>
          <Constant key="Parameter_4081" name="Ki1AKG" value="0.0150083"/>
          <Constant key="Parameter_4080" name="Ki1NADH" value="0.000330313"/>
          <Constant key="Parameter_4079" name="Ki2AKG" value="0.256252"/>
          <Constant key="Parameter_4078" name="Ki2NADH" value="0.0504"/>
          <Constant key="Parameter_4077" name="KiATP" value="0.579939"/>
          <Constant key="Parameter_4076" name="KmACCOA0" value="0.119974"/>
          <Constant key="Parameter_4075" name="KmOAA0" value="0.00123458"/>
          <Constant key="Parameter_4074" name="KmcsCIT" value="1.15715"/>
          <Constant key="Parameter_4073" name="KmcsCOA" value="9.6344e-005"/>
          <Constant key="Parameter_4072" name="Vmax" value="57.0584"/>
        </ListOfConstants>
        <KineticLaw function="Function_65">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_614">
              <SourceParameter reference="Metabolite_1"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_664">
              <SourceParameter reference="Metabolite_7"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_676">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_591">
              <SourceParameter reference="Metabolite_11"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_640">
              <SourceParameter reference="Metabolite_109"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_656">
              <SourceParameter reference="Parameter_4086"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_688">
              <SourceParameter reference="Parameter_4085"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_684">
              <SourceParameter reference="Parameter_4084"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_690">
              <SourceParameter reference="Parameter_4083"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_648">
              <SourceParameter reference="Parameter_4082"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_666">
              <SourceParameter reference="Parameter_4081"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_670">
              <SourceParameter reference="Parameter_4080"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_672">
              <SourceParameter reference="Parameter_4079"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_674">
              <SourceParameter reference="Parameter_4078"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_678">
              <SourceParameter reference="Parameter_4077"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_680">
              <SourceParameter reference="Parameter_4076"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_634">
              <SourceParameter reference="Parameter_4075"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_686">
              <SourceParameter reference="Parameter_4074"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_682">
              <SourceParameter reference="Parameter_4073"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_668">
              <SourceParameter reference="Metabolite_39"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_519">
              <SourceParameter reference="Metabolite_45"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_630">
              <SourceParameter reference="Parameter_4072"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_26" name="ACN_1" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_26">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T20:01:48Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_11" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_3" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_31" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4071" name="Keq" value="0.384503"/>
          <Constant key="Parameter_4070" name="KmACO" value="0.02001"/>
          <Constant key="Parameter_4069" name="KmCIT" value="0.0628882"/>
          <Constant key="Parameter_4068" name="KmICIT" value="9.31352"/>
          <Constant key="Parameter_4067" name="Vmax" value="9.72413"/>
        </ListOfConstants>
        <KineticLaw function="Function_66">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_683">
              <SourceParameter reference="Metabolite_3"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_687">
              <SourceParameter reference="Metabolite_11"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_667">
              <SourceParameter reference="Metabolite_31"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_679">
              <SourceParameter reference="Parameter_4071"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_671">
              <SourceParameter reference="ModelValue_6"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_675">
              <SourceParameter reference="ModelValue_5"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_663">
              <SourceParameter reference="ModelValue_4"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_691">
              <SourceParameter reference="Parameter_4067"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_27" name="ACN_2" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_27">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:56:27Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_3" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_31" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_11" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4066" name="Keq" value="3.49183"/>
          <Constant key="Parameter_4065" name="KmACO" value="0.02001"/>
          <Constant key="Parameter_4064" name="KmCIT" value="0.0628882"/>
          <Constant key="Parameter_4063" name="KmICIT" value="9.31352"/>
          <Constant key="Parameter_4062" name="Vmax" value="9.86571"/>
        </ListOfConstants>
        <KineticLaw function="Function_67">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_669">
              <SourceParameter reference="Metabolite_3"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_606">
              <SourceParameter reference="Metabolite_11"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_677">
              <SourceParameter reference="Metabolite_31"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_685">
              <SourceParameter reference="Parameter_4066"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_652">
              <SourceParameter reference="ModelValue_6"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_692">
              <SourceParameter reference="ModelValue_5"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_559">
              <SourceParameter reference="ModelValue_4"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_660">
              <SourceParameter reference="Parameter_4062"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_28" name="ICD" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_28">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T09:23:45Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_31" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_41" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_7" stoichiometry="1"/>
          <Product metabolite="Metabolite_43" stoichiometry="1"/>
          <Product metabolite="Metabolite_111" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_89" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4061" name="Keq" value="28.1631"/>
          <Constant key="Parameter_4060" name="KmAKG" value="0.038038"/>
          <Constant key="Parameter_4059" name="KmICIT" value="0.010989"/>
          <Constant key="Parameter_4058" name="KmNADP" value="0.005994"/>
          <Constant key="Parameter_4057" name="KmNADPH" value="0.000683333"/>
          <Constant key="Parameter_4056" name="kcat" value="2461.97"/>
        </ListOfConstants>
        <KineticLaw function="Function_68">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_694">
              <SourceParameter reference="Metabolite_7"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_689">
              <SourceParameter reference="Metabolite_31"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_698">
              <SourceParameter reference="Parameter_4061"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_704">
              <SourceParameter reference="Parameter_4060"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_700">
              <SourceParameter reference="Parameter_4059"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_702">
              <SourceParameter reference="Parameter_4058"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_706">
              <SourceParameter reference="Parameter_4057"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_673">
              <SourceParameter reference="Metabolite_41"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_696">
              <SourceParameter reference="Metabolite_43"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_693">
              <SourceParameter reference="Metabolite_89"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_487">
              <SourceParameter reference="Parameter_4056"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_29" name="LPD" reversible="false" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_29">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T17:04:37Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_109" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_7" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_37" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_39" stoichiometry="1"/>
          <Product metabolite="Metabolite_69" stoichiometry="1"/>
          <Product metabolite="Metabolite_111" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4055" name="KdAKG" value="14.9386"/>
          <Constant key="Parameter_4054" name="KmAKG" value="0.020002"/>
          <Constant key="Parameter_4053" name="KmCOA" value="0.0760076"/>
          <Constant key="Parameter_4052" name="KmNAD" value="0.0980098"/>
          <Constant key="Parameter_4051" name="Vmax" value="0.0684413"/>
          <Constant key="Parameter_4050" name="alpha" value="16.4304"/>
        </ListOfConstants>
        <KineticLaw function="Function_69">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_699">
              <SourceParameter reference="Metabolite_7"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_703">
              <SourceParameter reference="Metabolite_109"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_681">
              <SourceParameter reference="Parameter_4055"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_709">
              <SourceParameter reference="Parameter_4054"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_622">
              <SourceParameter reference="Parameter_4053"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_711">
              <SourceParameter reference="Parameter_4052"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_695">
              <SourceParameter reference="Metabolite_37"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_707">
              <SourceParameter reference="Parameter_4051"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_713">
              <SourceParameter reference="Parameter_4050"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_30" name="SK" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_30">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T16:20:56Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_101" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_69" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_113" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_105" stoichiometry="1"/>
          <Product metabolite="Metabolite_109" stoichiometry="1"/>
          <Product metabolite="Metabolite_67" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4049" name="Keq" value="1.15994"/>
          <Constant key="Parameter_4048" name="KmADP" value="0.00868475"/>
          <Constant key="Parameter_4047" name="KmATP" value="0.102321"/>
          <Constant key="Parameter_4046" name="KmCOA" value="0.255019"/>
          <Constant key="Parameter_4045" name="KmP" value="0.914709"/>
          <Constant key="Parameter_4044" name="KmSUC" value="0.800744"/>
          <Constant key="Parameter_4043" name="KmSUCCOA" value="0.0085"/>
          <Constant key="Parameter_4042" name="Vmax" value="76.8163"/>
        </ListOfConstants>
        <KineticLaw function="Function_70">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_710">
              <SourceParameter reference="Metabolite_101"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_705">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_716">
              <SourceParameter reference="Metabolite_109"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_720">
              <SourceParameter reference="Parameter_4049"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_722">
              <SourceParameter reference="Parameter_4048"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_728">
              <SourceParameter reference="Parameter_4047"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_730">
              <SourceParameter reference="Parameter_4046"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_726">
              <SourceParameter reference="Parameter_4045"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_732">
              <SourceParameter reference="Parameter_4044"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_724">
              <SourceParameter reference="Parameter_4043"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_697">
              <SourceParameter reference="Metabolite_113"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_718">
              <SourceParameter reference="Metabolite_67"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_644">
              <SourceParameter reference="Metabolite_69"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_714">
              <SourceParameter reference="Parameter_4042"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_31" name="SDH" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_31">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T15:38:25Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_135" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_67" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_21" stoichiometry="1"/>
          <Product metabolite="Metabolite_137" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4041" name="KefFUM" value="0.067048"/>
          <Constant key="Parameter_4040" name="KefSUC" value="0.0322425"/>
          <Constant key="Parameter_4039" name="Keq" value="2250"/>
          <Constant key="Parameter_4038" name="KmFUM" value="1.36019"/>
          <Constant key="Parameter_4037" name="KmQ" value="0.00160718"/>
          <Constant key="Parameter_4036" name="KmQH2" value="0.00611663"/>
          <Constant key="Parameter_4035" name="KmSUC" value="0.805727"/>
          <Constant key="Parameter_4034" name="Vmax" value="1.56184"/>
        </ListOfConstants>
        <KineticLaw function="Function_71">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_721">
              <SourceParameter reference="Metabolite_21"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_736">
              <SourceParameter reference="Parameter_4041"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_708">
              <SourceParameter reference="Parameter_4040"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_701">
              <SourceParameter reference="Parameter_4039"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_740">
              <SourceParameter reference="Parameter_4038"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_734">
              <SourceParameter reference="Parameter_4037"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_742">
              <SourceParameter reference="Parameter_4036"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_738">
              <SourceParameter reference="Parameter_4035"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_725">
              <SourceParameter reference="Metabolite_135"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_717">
              <SourceParameter reference="Metabolite_137"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_729">
              <SourceParameter reference="Metabolite_67"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_733">
              <SourceParameter reference="Parameter_4034"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_32" name="FUMA" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_32">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-07T12:23:44Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_21" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_35" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4033" name="Keq" value="4.99966"/>
          <Constant key="Parameter_4032" name="KmFUM" value="0.6"/>
          <Constant key="Parameter_4031" name="KmMAL" value="0.7"/>
          <Constant key="Parameter_4030" name="Vmax" value="53.3414"/>
        </ListOfConstants>
        <KineticLaw function="Function_72">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_739">
              <SourceParameter reference="Metabolite_21"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_665">
              <SourceParameter reference="Parameter_4033"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_719">
              <SourceParameter reference="Parameter_4032"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_727">
              <SourceParameter reference="Parameter_4031"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_735">
              <SourceParameter reference="Metabolite_35"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_743">
              <SourceParameter reference="Parameter_4030"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_33" name="MQO" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_33">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T16:57:52Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_35" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_57" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_45" stoichiometry="1"/>
          <Product metabolite="Metabolite_59" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4029" name="Keq" value="9"/>
          <Constant key="Parameter_4028" name="KmMAL" value="0.435"/>
          <Constant key="Parameter_4027" name="KmOAA" value="75.8036"/>
          <Constant key="Parameter_4026" name="KmQ" value="0.0414"/>
          <Constant key="Parameter_4025" name="KmQH2" value="8.77942"/>
          <Constant key="Parameter_4024" name="Vmax" value="4.62283"/>
        </ListOfConstants>
        <KineticLaw function="Function_73">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_748">
              <SourceParameter reference="Parameter_4029"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_750">
              <SourceParameter reference="Parameter_4028"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_754">
              <SourceParameter reference="Parameter_4027"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_752">
              <SourceParameter reference="Parameter_4026"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_756">
              <SourceParameter reference="Parameter_4025"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_715">
              <SourceParameter reference="Metabolite_35"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_744">
              <SourceParameter reference="Metabolite_45"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_737">
              <SourceParameter reference="Metabolite_57"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_746">
              <SourceParameter reference="Metabolite_59"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_731">
              <SourceParameter reference="Parameter_4024"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_34" name="MDH" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_34">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T16:59:32Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_59" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_45" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_35" stoichiometry="1"/>
          <Product metabolite="Metabolite_57" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_37" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_39" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4023" name="Keq" value="100000"/>
          <Constant key="Parameter_4022" name="KiNAD" value="0.0233122"/>
          <Constant key="Parameter_4021" name="KiNADH" value="0.000196981"/>
          <Constant key="Parameter_4020" name="KiOAA" value="2.46446"/>
          <Constant key="Parameter_4019" name="KmMAL" value="0.86"/>
          <Constant key="Parameter_4018" name="KmNAD" value="0.64"/>
          <Constant key="Parameter_4017" name="KmNADH" value="0.003"/>
          <Constant key="Parameter_4016" name="KmOAA" value="0.001"/>
          <Constant key="Parameter_4015" name="Vmax" value="6.11492"/>
        </ListOfConstants>
        <KineticLaw function="Function_74">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_758">
              <SourceParameter reference="Parameter_4023"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_766">
              <SourceParameter reference="Parameter_4022"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_760">
              <SourceParameter reference="Parameter_4021"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_772">
              <SourceParameter reference="Parameter_4020"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_768">
              <SourceParameter reference="Parameter_4019"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_764">
              <SourceParameter reference="Parameter_4018"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_770">
              <SourceParameter reference="Parameter_4017"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_762">
              <SourceParameter reference="Parameter_4016"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_745">
              <SourceParameter reference="Metabolite_35"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_712">
              <SourceParameter reference="Metabolite_37"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_753">
              <SourceParameter reference="Metabolite_39"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_749">
              <SourceParameter reference="Metabolite_45"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_757">
              <SourceParameter reference="Parameter_4015"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_35" name="ACEA" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_35">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T10:03:15Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_31" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_29" stoichiometry="1"/>
          <Product metabolite="Metabolite_67" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_47" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_51" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4014" name="KdICITsuc" value="0.00489074"/>
          <Constant key="Parameter_4013" name="KdPEP" value="1.05105"/>
          <Constant key="Parameter_4012" name="KdPEPglx" value="0.0358555"/>
          <Constant key="Parameter_4011" name="KdPEPicit" value="0.164263"/>
          <Constant key="Parameter_4010" name="KdPGA3" value="0.8"/>
          <Constant key="Parameter_4009" name="KdSUC" value="0.53053"/>
          <Constant key="Parameter_4008" name="Keq" value="8.8088"/>
          <Constant key="Parameter_4007" name="KmGLX" value="0.13013"/>
          <Constant key="Parameter_4006" name="KmICIT" value="0.063063"/>
          <Constant key="Parameter_4005" name="KmSUC" value="0.58941"/>
          <Constant key="Parameter_4004" name="Vmax" value="1.52595"/>
        </ListOfConstants>
        <KineticLaw function="Function_75">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_765">
              <SourceParameter reference="Metabolite_29"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_769">
              <SourceParameter reference="Metabolite_31"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_779">
              <SourceParameter reference="Parameter_4014"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_787">
              <SourceParameter reference="Parameter_4013"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_785">
              <SourceParameter reference="Parameter_4012"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_775">
              <SourceParameter reference="Parameter_4011"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_791">
              <SourceParameter reference="Parameter_4010"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_777">
              <SourceParameter reference="Parameter_4009"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_723">
              <SourceParameter reference="Parameter_4008"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_783">
              <SourceParameter reference="Parameter_4007"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_747">
              <SourceParameter reference="Parameter_4006"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_781">
              <SourceParameter reference="Parameter_4005"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_755">
              <SourceParameter reference="Metabolite_47"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_789">
              <SourceParameter reference="Metabolite_51"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_761">
              <SourceParameter reference="Metabolite_67"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_773">
              <SourceParameter reference="Parameter_4004"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_36" name="ACEB" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_36">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T11:32:16Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_1" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_29" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_109" stoichiometry="1"/>
          <Product metabolite="Metabolite_35" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4003" name="Keq" value="230000"/>
          <Constant key="Parameter_4002" name="KmACCOA" value="0.009"/>
          <Constant key="Parameter_4001" name="KmCOA" value="10.5652"/>
          <Constant key="Parameter_4000" name="KmGLX" value="0.021"/>
          <Constant key="Parameter_3999" name="KmMAL" value="15.0572"/>
          <Constant key="Parameter_3998" name="Vmax" value="0.352769"/>
        </ListOfConstants>
        <KineticLaw function="Function_76">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_788">
              <SourceParameter reference="Metabolite_1"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_780">
              <SourceParameter reference="Metabolite_109"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_784">
              <SourceParameter reference="Metabolite_29"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_751">
              <SourceParameter reference="Parameter_4003"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_759">
              <SourceParameter reference="Parameter_4002"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_793">
              <SourceParameter reference="Parameter_4001"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_767">
              <SourceParameter reference="Parameter_4000"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_795">
              <SourceParameter reference="Parameter_3999"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_776">
              <SourceParameter reference="Metabolite_35"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_792">
              <SourceParameter reference="Parameter_3998"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_37" name="ACEK_1" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_37">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:48:01Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_105" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_89" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_101" stoichiometry="1"/>
          <Product metabolite="Metabolite_91" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3997" name="Keq" value="888"/>
          <Constant key="Parameter_3996" name="k" value="1.25457"/>
        </ListOfConstants>
        <KineticLaw function="Function_77">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_778">
              <SourceParameter reference="Metabolite_101"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_771">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_797">
              <SourceParameter reference="Parameter_3997"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_741">
              <SourceParameter reference="Metabolite_89"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_786">
              <SourceParameter reference="Metabolite_91"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_796">
              <SourceParameter reference="Parameter_3996"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_38" name="ACEK_2" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_38">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:48:23Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_91" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_89" stoichiometry="1"/>
          <Product metabolite="Metabolite_113" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3995" name="Keq" value="400"/>
          <Constant key="Parameter_3994" name="k" value="0.0332"/>
        </ListOfConstants>
        <KineticLaw function="Function_78">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_801">
              <SourceParameter reference="Parameter_3995"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_799">
              <SourceParameter reference="Metabolite_113"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_763">
              <SourceParameter reference="Metabolite_89"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_782">
              <SourceParameter reference="Metabolite_91"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_798">
              <SourceParameter reference="Parameter_3994"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_39" name="EDD" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_39">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T19:07:17Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_53" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_33" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3993" name="Keq" value="1000"/>
          <Constant key="Parameter_3992" name="KmKDPG" value="0.318316"/>
          <Constant key="Parameter_3991" name="KmPGN" value="0.6"/>
          <Constant key="Parameter_3990" name="Vmax" value="0.111359"/>
        </ListOfConstants>
        <KineticLaw function="Function_79">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_790">
              <SourceParameter reference="Metabolite_33"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_804">
              <SourceParameter reference="Parameter_3993"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_808">
              <SourceParameter reference="Parameter_3992"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_806">
              <SourceParameter reference="Parameter_3991"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_794">
              <SourceParameter reference="Metabolite_53"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_802">
              <SourceParameter reference="Parameter_3990"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_40" name="EDA" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_40">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T19:06:34Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_33" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_25" stoichiometry="1"/>
          <Product metabolite="Metabolite_55" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3989" name="Keq" value="0.5"/>
          <Constant key="Parameter_3988" name="KmGAP" value="86.678"/>
          <Constant key="Parameter_3987" name="KmKDPG" value="0.06"/>
          <Constant key="Parameter_3986" name="KmPYR" value="10"/>
          <Constant key="Parameter_3985" name="Vmax" value="0.0775241"/>
        </ListOfConstants>
        <KineticLaw function="Function_80">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_774">
              <SourceParameter reference="Metabolite_25"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_805">
              <SourceParameter reference="Metabolite_33"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_812">
              <SourceParameter reference="Parameter_3989"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_816">
              <SourceParameter reference="Parameter_3988"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_814">
              <SourceParameter reference="Parameter_3987"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_818">
              <SourceParameter reference="Parameter_3986"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_810">
              <SourceParameter reference="Metabolite_55"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_809">
              <SourceParameter reference="Parameter_3985"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_41" name="NDHI" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_41">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T16:30:22Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_39" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_57" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_129" stoichiometry="4"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_37" stoichiometry="1"/>
          <Product metabolite="Metabolite_59" stoichiometry="1"/>
          <Product metabolite="Metabolite_147" stoichiometry="4"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_129" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_147" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_3984" name="Vmax" value="23.0735"/>
          <Constant key="Parameter_3983" name="Keq" value="27.6193"/>
        </ListOfConstants>
        <KineticLaw function="Function_93">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_903">
              <SourceParameter reference="Parameter_3984"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_895">
              <SourceParameter reference="Metabolite_129"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_899">
              <SourceParameter reference="Metabolite_147"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_892">
              <SourceParameter reference="Metabolite_39"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_905">
              <SourceParameter reference="Metabolite_57"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_907">
              <SourceParameter reference="Metabolite_37"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_909">
              <SourceParameter reference="Metabolite_59"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_911">
              <SourceParameter reference="ModelValue_7"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_42" name="PNT" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_42">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:54:17Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_37" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_43" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_39" stoichiometry="1"/>
          <Product metabolite="Metabolite_41" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3982" name="Keq" value="0.181891"/>
          <Constant key="Parameter_3981" name="k" value="2.49441"/>
        </ListOfConstants>
        <KineticLaw function="Function_81">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_822">
              <SourceParameter reference="Parameter_3982"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_815">
              <SourceParameter reference="Metabolite_37"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_803">
              <SourceParameter reference="Metabolite_39"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_820">
              <SourceParameter reference="Metabolite_41"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_811">
              <SourceParameter reference="Metabolite_43"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_819">
              <SourceParameter reference="Parameter_3981"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_43" name="ADK" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_43">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:48:35Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_103" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_105" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_101" stoichiometry="2"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3980" name="Keq" value="0.962758"/>
          <Constant key="Parameter_3979" name="k" value="0.242256"/>
        </ListOfConstants>
        <KineticLaw function="Function_82">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_824">
              <SourceParameter reference="Metabolite_101"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_807">
              <SourceParameter reference="Metabolite_103"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_813">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_826">
              <SourceParameter reference="Parameter_3980"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_823">
              <SourceParameter reference="Parameter_3979"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_44" name="ATP_SYN" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_44">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:48:42Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_101" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_113" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_147" stoichiometry="4"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_105" stoichiometry="1"/>
          <Product metabolite="Metabolite_129" stoichiometry="4"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_129" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_147" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_3978" name="Vmax" value="108.733"/>
          <Constant key="Parameter_3977" name="Keq" value="49.8315"/>
        </ListOfConstants>
        <KineticLaw function="Function_92">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_894">
              <SourceParameter reference="Parameter_3978"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_884">
              <SourceParameter reference="Metabolite_129"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_878">
              <SourceParameter reference="Metabolite_147"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_896">
              <SourceParameter reference="Metabolite_101"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_898">
              <SourceParameter reference="Metabolite_113"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_900">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_902">
              <SourceParameter reference="Parameter_3977"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_45" name="CYA" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_45">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:49:01Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_105" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_107" stoichiometry="1"/>
          <Product metabolite="Metabolite_113" stoichiometry="2"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_79" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_3976" name="Keq" value="2591.19"/>
          <Constant key="Parameter_3975" name="k" value="0.00414418"/>
          <Constant key="Parameter_3974" name="KaeiiaP" value="0.180981"/>
        </ListOfConstants>
        <KineticLaw function="Function_83">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_817">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_821">
              <SourceParameter reference="Metabolite_107"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_831">
              <SourceParameter reference="Parameter_3976"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_829">
              <SourceParameter reference="Metabolite_113"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_827">
              <SourceParameter reference="Parameter_3975"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_833">
              <SourceParameter reference="Metabolite_79"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_835">
              <SourceParameter reference="Parameter_3974"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_46" name="DOS" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_46">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:49:37Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_107" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_103" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3973" name="k" value="0.00828"/>
          <Constant key="Parameter_3972" name="Keq" value="674.242"/>
        </ListOfConstants>
        <KineticLaw function="Function_84">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_828">
              <SourceParameter reference="Metabolite_103"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_832">
              <SourceParameter reference="Metabolite_107"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_825">
              <SourceParameter reference="Parameter_3972"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_836">
              <SourceParameter reference="Parameter_3973"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_47" name="ACK" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_47">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:48:31Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_5" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_101" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_139" stoichiometry="1"/>
          <Product metabolite="Metabolite_105" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3971" name="Keq" value="174"/>
          <Constant key="Parameter_3970" name="KmACE" value="7"/>
          <Constant key="Parameter_3969" name="KmACP" value="0.16"/>
          <Constant key="Parameter_3968" name="KmADP" value="0.5"/>
          <Constant key="Parameter_3967" name="KmATP" value="0.07"/>
          <Constant key="Parameter_3966" name="Vmax" value="7.23"/>
        </ListOfConstants>
        <KineticLaw function="Function_85">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_842">
              <SourceParameter reference="Metabolite_139"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_838">
              <SourceParameter reference="Metabolite_5"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_840">
              <SourceParameter reference="Metabolite_101"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_844">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_846">
              <SourceParameter reference="Parameter_3971"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_852">
              <SourceParameter reference="Parameter_3970"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_848">
              <SourceParameter reference="Parameter_3969"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_850">
              <SourceParameter reference="Parameter_3968"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_854">
              <SourceParameter reference="Parameter_3967"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_830">
              <SourceParameter reference="Parameter_3966"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_837">
              <SourceParameter reference="Compartment_1"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_48" name="ACS" reversible="false" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_48">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:59:08Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_139" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_105" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_109" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_1" stoichiometry="1"/>
          <Product metabolite="Metabolite_103" stoichiometry="1"/>
          <Product metabolite="Metabolite_113" stoichiometry="2"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3965" name="KmACE" value="0.07"/>
          <Constant key="Parameter_3964" name="KmATP" value="0.1"/>
          <Constant key="Parameter_3963" name="KmCOA" value="0.01"/>
          <Constant key="Parameter_3962" name="Vmax" value="7.3"/>
        </ListOfConstants>
        <KineticLaw function="Function_86">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_847">
              <SourceParameter reference="Metabolite_139"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_843">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_839">
              <SourceParameter reference="Metabolite_109"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_800">
              <SourceParameter reference="Parameter_3965"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_857">
              <SourceParameter reference="Parameter_3964"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_859">
              <SourceParameter reference="Parameter_3963"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_851">
              <SourceParameter reference="Parameter_3962"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_855">
              <SourceParameter reference="Compartment_1"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_49" name="PTA" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_49">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T17:09:21Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_1" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_113" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_109" stoichiometry="1"/>
          <Product metabolite="Metabolite_5" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3961" name="Keq" value="0.0053952"/>
          <Constant key="Parameter_3960" name="KiACCOA" value="0.2"/>
          <Constant key="Parameter_3959" name="KiACP" value="0.2"/>
          <Constant key="Parameter_3958" name="KiCOA" value="0.029"/>
          <Constant key="Parameter_3957" name="KiP" value="13.54"/>
          <Constant key="Parameter_3956" name="KmACP" value="0.7"/>
          <Constant key="Parameter_3955" name="KmP" value="6.1"/>
          <Constant key="Parameter_3954" name="Vmax" value="2.7"/>
        </ListOfConstants>
        <KineticLaw function="Function_87">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_856">
              <SourceParameter reference="Metabolite_1"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_849">
              <SourceParameter reference="Metabolite_5"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_861">
              <SourceParameter reference="Metabolite_109"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_863">
              <SourceParameter reference="Parameter_3961"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_865">
              <SourceParameter reference="Parameter_3960"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_871">
              <SourceParameter reference="Parameter_3959"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_873">
              <SourceParameter reference="Parameter_3958"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_869">
              <SourceParameter reference="Parameter_3957"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_875">
              <SourceParameter reference="Parameter_3956"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_867">
              <SourceParameter reference="Parameter_3955"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_841">
              <SourceParameter reference="Metabolite_113"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_860">
              <SourceParameter reference="Parameter_3954"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_50" name="PTS_0" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_50">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T12:14:31Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_73" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_47" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_75" stoichiometry="1"/>
          <Product metabolite="Metabolite_55" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3953" name="KmPEP" value="0.6"/>
          <Constant key="Parameter_3952" name="KmPYR" value="1"/>
          <Constant key="Parameter_3951" name="kF" value="12000"/>
          <Constant key="Parameter_3950" name="kR" value="8000"/>
        </ListOfConstants>
        <KineticLaw function="Function_88">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_864">
              <SourceParameter reference="Parameter_3953"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_879">
              <SourceParameter reference="Parameter_3952"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_868">
              <SourceParameter reference="Metabolite_47"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_877">
              <SourceParameter reference="Metabolite_55"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_872">
              <SourceParameter reference="Metabolite_73"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_834">
              <SourceParameter reference="Metabolite_75"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_876">
              <SourceParameter reference="Parameter_3951"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_853">
              <SourceParameter reference="Parameter_3950"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_51" name="PTS_1" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_51">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-03T00:49:35Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_85" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_75" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_87" stoichiometry="1"/>
          <Product metabolite="Metabolite_73" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3949" name="k1" value="200000"/>
          <Constant key="Parameter_3948" name="k2" value="8000"/>
        </ListOfConstants>
        <KineticLaw function="Function_14">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_62">
              <SourceParameter reference="Parameter_3949"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_39">
              <SourceParameter reference="Metabolite_85"/>
              <SourceParameter reference="Metabolite_75"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_67">
              <SourceParameter reference="Parameter_3948"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_71">
              <SourceParameter reference="Metabolite_87"/>
              <SourceParameter reference="Metabolite_73"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_52" name="PTS_2" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_52">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T10:24:08Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_77" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_87" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_79" stoichiometry="1"/>
          <Product metabolite="Metabolite_85" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3947" name="k1" value="61000"/>
          <Constant key="Parameter_3946" name="k2" value="47000"/>
        </ListOfConstants>
        <KineticLaw function="Function_14">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_62">
              <SourceParameter reference="Parameter_3947"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_39">
              <SourceParameter reference="Metabolite_77"/>
              <SourceParameter reference="Metabolite_87"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_67">
              <SourceParameter reference="Parameter_3946"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_71">
              <SourceParameter reference="Metabolite_79"/>
              <SourceParameter reference="Metabolite_85"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_53" name="PTS_3" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_53">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T12:14:41Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_81" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_79" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_83" stoichiometry="1"/>
          <Product metabolite="Metabolite_77" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3945" name="k1" value="11000"/>
          <Constant key="Parameter_3944" name="k2" value="4000"/>
        </ListOfConstants>
        <KineticLaw function="Function_14">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_62">
              <SourceParameter reference="Parameter_3945"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_39">
              <SourceParameter reference="Metabolite_81"/>
              <SourceParameter reference="Metabolite_79"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_67">
              <SourceParameter reference="Parameter_3944"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_71">
              <SourceParameter reference="Metabolite_83"/>
              <SourceParameter reference="Metabolite_77"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_54" name="PTS_4" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_54">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T12:14:34Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_149" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_83" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_23" stoichiometry="1"/>
          <Product metabolite="Metabolite_81" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3943" name="KmG6P" value="2125.91"/>
          <Constant key="Parameter_3942" name="KmGLC" value="0.02"/>
          <Constant key="Parameter_3941" name="kF" value="4000"/>
          <Constant key="Parameter_3940" name="kR" value="1e-005"/>
        </ListOfConstants>
        <KineticLaw function="Function_89">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_887">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_870">
              <SourceParameter reference="Metabolite_149"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_889">
              <SourceParameter reference="Parameter_3943"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_881">
              <SourceParameter reference="Parameter_3942"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_880">
              <SourceParameter reference="Compartment_1"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_885">
              <SourceParameter reference="Metabolite_81"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_862">
              <SourceParameter reference="Metabolite_83"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_858">
              <SourceParameter reference="Parameter_3941"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_883">
              <SourceParameter reference="Parameter_3940"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_55" name="_GLC_FEED" reversible="false" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_55">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:49:02Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfProducts>
          <Product metabolite="Metabolite_141" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3939" name="FEED" value="0.23"/>
        </ListOfConstants>
        <KineticLaw function="Function_90">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_886">
              <SourceParameter reference="ModelValue_0"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_890">
              <SourceParameter reference="Compartment_1"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_882">
              <SourceParameter reference="Compartment_3"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_56" name="CYTBO" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_56">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T09:53:39Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_59" stoichiometry="2"/>
          <Substrate metabolite="Metabolite_129" stoichiometry="8"/>
          <Substrate metabolite="Metabolite_133" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_57" stoichiometry="2"/>
          <Product metabolite="Metabolite_147" stoichiometry="8"/>
          <Product metabolite="Metabolite_131" stoichiometry="2"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_129" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_147" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_3938" name="Vmax" value="8.54045"/>
          <Constant key="Parameter_3937" name="Keq" value="12.067"/>
        </ListOfConstants>
        <KineticLaw function="Function_94">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_912">
              <SourceParameter reference="Parameter_3938"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_904">
              <SourceParameter reference="Metabolite_129"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_908">
              <SourceParameter reference="Metabolite_147"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_917">
              <SourceParameter reference="Parameter_3937"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_913">
              <SourceParameter reference="Metabolite_133"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_897">
              <SourceParameter reference="Metabolite_59"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_915">
              <SourceParameter reference="Metabolite_57"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_57" name="SQR" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_57">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T13:54:45Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_137" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_57" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_135" stoichiometry="1"/>
          <Product metabolite="Metabolite_59" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3936" name="Keq" value="0.94033"/>
          <Constant key="Parameter_3935" name="Vmax" value="3.41617"/>
        </ListOfConstants>
        <KineticLaw function="Function_95">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_914">
              <SourceParameter reference="Metabolite_137"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_866">
              <SourceParameter reference="Metabolite_57"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_910">
              <SourceParameter reference="Metabolite_135"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_920">
              <SourceParameter reference="Metabolite_59"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_922">
              <SourceParameter reference="Parameter_3936"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_918">
              <SourceParameter reference="Parameter_3935"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_58" name="NDHII" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_58">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-17T08:01:41Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_39" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_57" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_37" stoichiometry="1"/>
          <Product metabolite="Metabolite_59" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3934" name="Vmax" value="30.8306"/>
          <Constant key="Parameter_3933" name="Keq" value="27.6193"/>
        </ListOfConstants>
        <KineticLaw function="Function_96">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_923">
              <SourceParameter reference="Parameter_3934"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_919">
              <SourceParameter reference="Metabolite_39"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_901">
              <SourceParameter reference="Metabolite_57"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_924">
              <SourceParameter reference="Metabolite_37"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_926">
              <SourceParameter reference="Metabolite_59"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_928">
              <SourceParameter reference="ModelValue_7"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_59" name="GROWTH" reversible="false" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_59">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-21T14:35:01Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_23" stoichiometry="116"/>
          <Substrate metabolite="Metabolite_15" stoichiometry="204"/>
          <Substrate metabolite="Metabolite_51" stoichiometry="845"/>
          <Substrate metabolite="Metabolite_45" stoichiometry="1010"/>
          <Substrate metabolite="Metabolite_7" stoichiometry="610"/>
          <Substrate metabolite="Metabolite_55" stoichiometry="1601"/>
          <Substrate metabolite="Metabolite_61" stoichiometry="507"/>
          <Substrate metabolite="Metabolite_47" stoichiometry="293"/>
          <Substrate metabolite="Metabolite_25" stoichiometry="73"/>
          <Substrate metabolite="Metabolite_17" stoichiometry="40"/>
          <Substrate metabolite="Metabolite_43" stoichiometry="10169"/>
          <Substrate metabolite="Metabolite_1" stoichiometry="2118"/>
          <Substrate metabolite="Metabolite_37" stoichiometry="2004"/>
          <Substrate metabolite="Metabolite_105" stoichiometry="30508"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_41" stoichiometry="10169"/>
          <Product metabolite="Metabolite_109" stoichiometry="2118"/>
          <Product metabolite="Metabolite_39" stoichiometry="2004"/>
          <Product metabolite="Metabolite_101" stoichiometry="30508"/>
          <Product metabolite="Metabolite_113" stoichiometry="30508"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3932" name="Vmax" value="9.74137"/>
          <Constant key="Parameter_3931" name="KmG6P" value="1.20911"/>
          <Constant key="Parameter_3930" name="KmE4P" value="1.63298"/>
          <Constant key="Parameter_3929" name="KmPGA3" value="0.0764619"/>
          <Constant key="Parameter_3928" name="KmOAA" value="0.0248351"/>
          <Constant key="Parameter_3927" name="KmAKG" value="5.11989"/>
          <Constant key="Parameter_3926" name="KmPYR" value="0.00463904"/>
          <Constant key="Parameter_3925" name="KmR5P" value="0.021234"/>
          <Constant key="Parameter_3924" name="KmPEP" value="0.457734"/>
          <Constant key="Parameter_3923" name="KmGAP" value="0.024854"/>
          <Constant key="Parameter_3922" name="KmF6P" value="0.366423"/>
          <Constant key="Parameter_3921" name="KmNADPH" value="3.59774"/>
          <Constant key="Parameter_3920" name="KmACCOA" value="0.0494404"/>
          <Constant key="Parameter_3919" name="KmNAD" value="2.82239"/>
          <Constant key="Parameter_3918" name="KmATP" value="0.0468266"/>
        </ListOfConstants>
        <KineticLaw function="Function_97">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_929">
              <SourceParameter reference="Parameter_3932"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_925">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_906">
              <SourceParameter reference="Metabolite_15"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_930">
              <SourceParameter reference="Metabolite_51"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_932">
              <SourceParameter reference="Metabolite_45"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_934">
              <SourceParameter reference="Metabolite_7"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_936">
              <SourceParameter reference="Metabolite_55"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_938">
              <SourceParameter reference="Metabolite_61"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_940">
              <SourceParameter reference="Metabolite_47"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_942">
              <SourceParameter reference="Metabolite_25"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_944">
              <SourceParameter reference="Metabolite_17"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_946">
              <SourceParameter reference="Metabolite_43"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_948">
              <SourceParameter reference="Metabolite_1"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_950">
              <SourceParameter reference="Metabolite_37"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_952">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_954">
              <SourceParameter reference="Parameter_3931"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_956">
              <SourceParameter reference="Parameter_3930"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_958">
              <SourceParameter reference="Parameter_3929"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_960">
              <SourceParameter reference="Parameter_3928"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_962">
              <SourceParameter reference="Parameter_3927"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_964">
              <SourceParameter reference="Parameter_3926"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_966">
              <SourceParameter reference="Parameter_3925"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_968">
              <SourceParameter reference="Parameter_3924"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_970">
              <SourceParameter reference="Parameter_3923"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_972">
              <SourceParameter reference="Parameter_3922"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_974">
              <SourceParameter reference="Parameter_3921"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_976">
              <SourceParameter reference="Parameter_3920"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_978">
              <SourceParameter reference="Parameter_3919"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_980">
              <SourceParameter reference="Parameter_3918"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_60" name="ATP_MAINTENANCE" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_60">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-21T14:43:37Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_105" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_101" stoichiometry="1"/>
          <Product metabolite="Metabolite_113" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3917" name="Vmax" value="1.30166"/>
          <Constant key="Parameter_3916" name="Keq" value="3.63369"/>
        </ListOfConstants>
        <KineticLaw function="Function_91">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_874">
              <SourceParameter reference="Parameter_3917"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_845">
              <SourceParameter reference="Metabolite_101"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_891">
              <SourceParameter reference="Metabolite_113"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_888">
              <SourceParameter reference="Metabolite_105"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_893">
              <SourceParameter reference="Parameter_3916"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_61" name="XCH_GLC" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_61">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-23T16:07:48Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_141" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_149" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3915" name="Vmax" value="100"/>
          <Constant key="Parameter_3914" name="Km" value="10"/>
        </ListOfConstants>
        <KineticLaw function="Function_100">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_979">
              <SourceParameter reference="Parameter_3915"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_963">
              <SourceParameter reference="Metabolite_141"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_945">
              <SourceParameter reference="Metabolite_149"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_949">
              <SourceParameter reference="Parameter_3914"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_62" name="PIT" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_62">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-23T17:08:39Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_151" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_147" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_113" stoichiometry="1"/>
          <Product metabolite="Metabolite_129" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_147" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_129" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_3913" name="Vmax" value="7.146"/>
          <Constant key="Parameter_3912" name="KmP" value="12.18"/>
          <Constant key="Parameter_3911" name="KmPp" value="0.025"/>
          <Constant key="Parameter_3910" name="Kr" value="0.1"/>
        </ListOfConstants>
        <KineticLaw function="Function_98">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_981">
              <SourceParameter reference="Parameter_3913"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_977">
              <SourceParameter reference="Metabolite_147"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_973">
              <SourceParameter reference="Metabolite_129"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_969">
              <SourceParameter reference="Metabolite_151"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_965">
              <SourceParameter reference="Parameter_3911"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_961">
              <SourceParameter reference="Parameter_3910"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_957">
              <SourceParameter reference="Metabolite_113"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_953">
              <SourceParameter reference="Parameter_3912"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_63" name="XCH_P" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_63">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-31T12:45:19Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_143" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_151" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3909" name="P" value="1.50545e+021"/>
          <Constant key="Parameter_3908" name="Vmax" value="100"/>
          <Constant key="Parameter_3907" name="Km" value="10"/>
        </ListOfConstants>
        <KineticLaw function="Function_100">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_979">
              <SourceParameter reference="Parameter_3908"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_963">
              <SourceParameter reference="Metabolite_143"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_945">
              <SourceParameter reference="Metabolite_151"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_949">
              <SourceParameter reference="Parameter_3907"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_64" name="XCH_ACE1" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_64">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-03T00:04:32Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_139" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_153" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3906" name="P" value="2.8108e+016"/>
          <Constant key="Parameter_3905" name="Vmax" value="100"/>
          <Constant key="Parameter_3904" name="Km" value="10"/>
        </ListOfConstants>
        <KineticLaw function="Function_100">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_979">
              <SourceParameter reference="Parameter_3905"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_963">
              <SourceParameter reference="Metabolite_139"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_945">
              <SourceParameter reference="Metabolite_153"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_949">
              <SourceParameter reference="Parameter_3904"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_65" name="_ACE_OUT" reversible="false" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_65">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-03T00:06:14Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_145" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfConstants>
          <Constant key="Parameter_3903" name="k1" value="5.556e-005"/>
        </ListOfConstants>
        <KineticLaw function="Function_13">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_81">
              <SourceParameter reference="Parameter_3903"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_79">
              <SourceParameter reference="Metabolite_145"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_66" name="XCH_ACE2" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_66">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-03T00:13:46Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_153" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_145" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3902" name="P" value="1.12426e+019"/>
          <Constant key="Parameter_3901" name="Vmax" value="100"/>
          <Constant key="Parameter_3900" name="Km" value="10"/>
        </ListOfConstants>
        <KineticLaw function="Function_100">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_979">
              <SourceParameter reference="Parameter_3901"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_963">
              <SourceParameter reference="Metabolite_153"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_945">
              <SourceParameter reference="Metabolite_145"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_949">
              <SourceParameter reference="Parameter_3900"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_67" name="GL6P_HYDROLYSIS" reversible="true" fast="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_67">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-22T11:03:40Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_27" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_53" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_3899" name="KGl6Phydrol" value="0.000167"/>
          <Constant key="Parameter_3898" name="KeqGl6Phydrol" value="42.8"/>
        </ListOfConstants>
        <KineticLaw function="Function_99">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_951">
              <SourceParameter reference="Parameter_3899"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_959">
              <SourceParameter reference="Metabolite_27"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_967">
              <SourceParameter reference="Metabolite_53"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_975">
              <SourceParameter reference="Parameter_3898"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
    </ListOfReactions>
    <ListOfModelParameterSets activeSet="ModelParameterSet_1">
      <ModelParameterSet key="ModelParameterSet_1" name="Initial State">
        <ModelParameterGroup cn="String=Initial Time" type="Group">
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism" value="0" type="Model" simulationType="time"/>
        </ModelParameterGroup>
        <ModelParameterGroup cn="String=Initial Compartment Sizes" type="Group">
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]" value="1" type="Compartment" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[extracellular]" value="100" type="Compartment" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_periplasm]" value="0.25" type="Compartment" simulationType="fixed"/>
        </ModelParameterGroup>
        <ModelParameterGroup cn="String=Initial Species Values" type="Group">
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ACCOA]" value="9.318843183080258e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ACO]" value="1.936106650655732e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ACP]" value="8.666420504823286e+018" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[AKG]" value="3.600459862224256e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[BPG]" value="3.939117017287109e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[CIT]" value="5.39165207053416e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[DAP]" value="2.632242451004475e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[E4P]" value="7.907068085616183e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[F6P]" value="1.576394850846697e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[FDP]" value="1.697088696281518e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[FUM]" value="1.283408617631331e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[G6P]" value="5.185843963253621e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[GAP]" value="7.05693301383649e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[GL6P]" value="1.964213101596648e+018" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[GLX]" value="5.305754012188268e+018" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ICIT]" value="6.370367694157086e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[KDPG]" value="5.228284334351698e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MAL]" value="6.215772458634581e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[NAD]" value="8.500516016790164e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[NADH]" value="9.542465977025343e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[NADP]" value="1.010739458841187e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[NADPH]" value="5.369509789454854e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[OAA]" value="7.698681357565698e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[PEP]" value="6.004306279627797e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[PGA2]" value="2.278158803291449e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[PGA3]" value="4.193061747813792e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[PGN]" value="7.925125419323508e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[PYR]" value="1.426593487227837e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[Q]" value="1.963925127694926e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[QH2]" value="4.058216660206235e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[R5P]" value="6.434157914056002e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[RU5P]" value="2.058528180838142e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[S7P]" value="8.550546593595767e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[SUC]" value="1.301115620269875e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[SUCCOA]" value="2.474364250168681e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[X5P]" value="3.047312090476589e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ei]" value="2.011473853979885e+017" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[eiP]" value="3.842826679553765e+018" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[eiia]" value="8.552558885904762e+018" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[eiiaP]" value="2.894051911487251e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[eiicb]" value="2.841906957912192e+016" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[eiicbP]" value="2.063314616344276e+017" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[hpr]" value="1.15150454342188e+017" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[hprP]" value="3.162113325697529e+018" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[icd]" value="3.268480204715989e+017" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[icdP]" value="3.098828960171077e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[tal]" value="1.674723810868397e+018" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[talC3]" value="3.445812690814797e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[tkt]" value="8.329093572866015e+018" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[tktC2]" value="3.382589890806846e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ADP]" value="3.603135976315973e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[AMP]" value="1.121639644055652e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ATP]" value="1.549016191114563e+021" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[CAMP]" value="5.55922809558007e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[COA]" value="3.011070895e+020" type="Species" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[HCO3]" value="8.430998506e+020" type="Species" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[P]" value="5.879987842207703e+021" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG]" value="6.022141789999999e+020" type="Species" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MgADP]" value="1.581911646485274e+020" type="Species" simulationType="assignment"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MgATP]" value="1.42797608318511e+021" type="Species" simulationType="assignment"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MgFDP]" value="2.492053885875944e+019" type="Species" simulationType="assignment"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ASP]" value="7.045905894300009e+020" type="Species" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[CYS]" value="5.118820521500001e+019" type="Species" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MN]" value="1.806642537000001e+020" type="Species" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[Hin]" value="1.904368444786942e+016" type="Species" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[H2O]" value="6.022141789999999e+020" type="Species" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[O2]" value="1.2646497759e+020" type="Species" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[FAD]" value="1.529107897748966e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[FADH2]" value="4.493033892713735e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ACE]" value="1.124380943115288e+017" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[extracellular],Vector=Metabolites[GLCx]" value="3.814045050430201e+020" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[extracellular],Vector=Metabolites[Px]" value="6.022141790000004e+023" type="Species" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[extracellular],Vector=Metabolites[ACEx]" value="1.124256015786816e+019" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_periplasm],Vector=Metabolites[Hout]" value="8091473706854326" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_periplasm],Vector=Metabolites[GLCp]" value="6.072381096825504e+017" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_periplasm],Vector=Metabolites[Pp]" value="1.505449345660231e+021" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_periplasm],Vector=Metabolites[ACEp]" value="2.810796198627631e+016" type="Species" simulationType="reactions"/>
        </ModelParameterGroup>
        <ModelParameterGroup cn="String=Initial Global Quantities" type="Group">
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[FEED]" value="0.23" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdADPMg]" value="1.27771" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdATPMg]" value="0.0847634" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdFDPMg]" value="5.81" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KmICIT_ACN]" value="9.313520000000001" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KmCIT_ACN]" value="0.06288820000000001" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KmACO_ACN]" value="0.02001" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KeqNDH]" value="27.6193" type="ModelValue" simulationType="fixed"/>
        </ModelParameterGroup>
        <ModelParameterGroup cn="String=Kinetic Parameters" type="Group">
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI],ParameterGroup=Parameters,Parameter=Keq" value="0.36" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI],ParameterGroup=Parameters,Parameter=KmF6P" value="0.147" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI],ParameterGroup=Parameters,Parameter=KmG6P" value="0.28" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI],ParameterGroup=Parameters,Parameter=KmPEP" value="1.999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI],ParameterGroup=Parameters,Parameter=Vmax" value="2.32456" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI],ParameterGroup=Parameters,Parameter=KmPGN" value="0.515958" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KefrADP" value="0.07352640000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KefrPEP" value="19.98" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KeftADP" value="9.009" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KeftPEP" value="0.26026" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=Keq" value="1998" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KirADP" value="54.945" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KirATP" value="2.4975e-005" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KirF6P" value="1.84615" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KirFDP" value="0.045954" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KitADP" value="80.08" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KitATP" value="0.014014" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KitF6P" value="0.008568559999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KitFDP" value="50.5505" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmrADP" value="0.690009" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmrATPMg" value="8.12187e-005" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmrF6P" value="2.05205e-005" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmrFDP" value="10.01" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmtADP" value="2.002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmtATPMg" value="3.34334" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmtF6P" value="32.967" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmtFDP" value="9.99" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=L0" value="14.0851" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=Vmax" value="0.185253" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=Wr" value="0.0237041" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=Wt" value="0.146735" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=n" value="4" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA],ParameterGroup=Parameters,Parameter=Keq" value="0.18981" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA],ParameterGroup=Parameters,Parameter=KmDAP" value="0.13001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA],ParameterGroup=Parameters,Parameter=KmFDP" value="0.12012" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA],ParameterGroup=Parameters,Parameter=KmGAP" value="0.13001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA],ParameterGroup=Parameters,Parameter=KmPEP" value="0.5" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA],ParameterGroup=Parameters,Parameter=Vmax" value="21.6978" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[TPI]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[TPI],ParameterGroup=Parameters,Parameter=Keq" value="0.270203" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[TPI],ParameterGroup=Parameters,Parameter=KmDAP" value="0.01" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[TPI],ParameterGroup=Parameters,Parameter=KmGAP" value="1.89301" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[TPI],ParameterGroup=Parameters,Parameter=Vmax" value="24.1843" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=Keq" value="20" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=KmBPG" value="0.2" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=KmGAP" value="2.47265" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=KmNAD" value="0.0110454" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=KmNADH" value="3.69797" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=KmP" value="0.017" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=Vmax" value="8.66573" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK],ParameterGroup=Parameters,Parameter=Keq" value="99.99250000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK],ParameterGroup=Parameters,Parameter=KmADPMg" value="0.08541600000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK],ParameterGroup=Parameters,Parameter=KmATPMg" value="3.47737" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK],ParameterGroup=Parameters,Parameter=KmBPG" value="0.0113296" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK],ParameterGroup=Parameters,Parameter=KmPGA3" value="2.45722" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK],ParameterGroup=Parameters,Parameter=Vmax" value="16.1089" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GPM]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GPM],ParameterGroup=Parameters,Parameter=Keq" value="0.565818" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GPM],ParameterGroup=Parameters,Parameter=KmPGA2" value="1.9153" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GPM],ParameterGroup=Parameters,Parameter=KmPGA3" value="0.115" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GPM],ParameterGroup=Parameters,Parameter=Vmax" value="10.9934" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ENO]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ENO],ParameterGroup=Parameters,Parameter=Keq" value="3" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ENO],ParameterGroup=Parameters,Parameter=KmPEP" value="0.1" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ENO],ParameterGroup=Parameters,Parameter=KmPGA2" value="0.1" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ENO],ParameterGroup=Parameters,Parameter=Vmax" value="11.7189" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrFDP" value="0.449149" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrG6P" value="0.158746" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrGL6P" value="0.150482" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrR5P" value="9.33254" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrRU5P" value="1.53591" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrS7P" value="0.0785955" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrX5P" value="0.677374" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KeftATP" value="3.69117" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KeftSUCCOA" value="8.264060000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KirADP" value="0.517585" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KirATP" value="96.0333" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KirPEP" value="0.181056" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KirPYR" value="15.1403" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KirPyrATP" value="230.781" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KitADP" value="0.224911" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KitATP" value="0.039564" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KitPEP" value="0.465672" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KitPYR" value="0.2499" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KitPyrATP" value="11.3691" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KmrADPMg" value="0.326144" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KmrPEP" value="5.563680000000001e-007" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KmtADPMg" value="0.054678" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KmtPEP" value="0.11475" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=L0" value="50.4818" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=Vmax" value="0.7471600000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=n" value="4" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ZWF]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ZWF],ParameterGroup=Parameters,Parameter=KdG6P" value="0.192" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ZWF],ParameterGroup=Parameters,Parameter=KdGL6P" value="0.02" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ZWF],ParameterGroup=Parameters,Parameter=Keq" value="60000000000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ZWF],ParameterGroup=Parameters,Parameter=KmG6P" value="0.118525" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ZWF],ParameterGroup=Parameters,Parameter=KmGL6P" value="0.328629" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ZWF],ParameterGroup=Parameters,Parameter=KmNADP" value="0.0274" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ZWF],ParameterGroup=Parameters,Parameter=KmNADPH" value="0.0168" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ZWF],ParameterGroup=Parameters,Parameter=Vmax" value="0.2658" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGL]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGL],ParameterGroup=Parameters,Parameter=Keq" value="42.7572" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGL],ParameterGroup=Parameters,Parameter=KiG6P" value="2.0001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGL],ParameterGroup=Parameters,Parameter=KmGL6P" value="0.022977" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGL],ParameterGroup=Parameters,Parameter=KmPGN" value="9.99" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGL],ParameterGroup=Parameters,Parameter=Vmax" value="11.5967" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KdHCO3" value="58.8951" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KdHCO3NADPH" value="9.7257" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KdNADP" value="0.116989" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KdNADPH" value="0.00340034" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KdRu5P" value="0.0440044" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KefATP" value="0.0650065" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KefFbP" value="0.0129987" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KefNADPATP" value="0.139986" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KefNADPFbP" value="0.00519948" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=Keq" value="49.962" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KmHCO3" value="6.41899" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KmNADP" value="0.049" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KmNADPH" value="68.3828" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KmPGN" value="0.093" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=KmRU5P" value="45.1977" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GND],ParameterGroup=Parameters,Parameter=Vmax" value="4.08105" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[RPE]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[RPE],ParameterGroup=Parameters,Parameter=Keq" value="1.5015" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[RPE],ParameterGroup=Parameters,Parameter=KmRU5P" value="0.872522" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[RPE],ParameterGroup=Parameters,Parameter=KmX5P" value="0.893607" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[RPE],ParameterGroup=Parameters,Parameter=Vmax" value="6.00103" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[RPI]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[RPI],ParameterGroup=Parameters,Parameter=Keq" value="0.330093" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[RPI],ParameterGroup=Parameters,Parameter=KmE4P" value="0.67067" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[RPI],ParameterGroup=Parameters,Parameter=KmR5P" value="3.09715" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[RPI],ParameterGroup=Parameters,Parameter=KmRU5P" value="4.40263" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[RPI],ParameterGroup=Parameters,Parameter=Vmax" value="8" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[X5P_GAP_TKT]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[X5P_GAP_TKT],ParameterGroup=Parameters,Parameter=Keq" value="1.001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[X5P_GAP_TKT],ParameterGroup=Parameters,Parameter=kcat" value="40" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[F6P_E4P_TKT]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[F6P_E4P_TKT],ParameterGroup=Parameters,Parameter=Keq" value="0.5004999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[F6P_E4P_TKT],ParameterGroup=Parameters,Parameter=kcat" value="40.0002" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[S7P_R5P_TKT]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[S7P_R5P_TKT],ParameterGroup=Parameters,Parameter=Keq" value="0.33033" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[S7P_R5P_TKT],ParameterGroup=Parameters,Parameter=kcat" value="199.979" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[F6P_GAP_TAL]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[F6P_GAP_TAL],ParameterGroup=Parameters,Parameter=Keq" value="0.11011" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[F6P_GAP_TAL],ParameterGroup=Parameters,Parameter=kcat" value="119.992" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[S7P_E4P_TAL]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[S7P_E4P_TAL],ParameterGroup=Parameters,Parameter=Keq" value="26.6266" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[S7P_E4P_TAL],ParameterGroup=Parameters,Parameter=kcat" value="99.9991" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirAMP" value="0.00122122" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirAMPFDP" value="0.256256" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirF6P" value="1.12112" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirF6PMg" value="0.384615" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirFDP" value="1.35327" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirFDPMg" value="0.75924" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirFDPMgMg" value="0.356356" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirP" value="3.16316" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirPF6P" value="6.60538" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirPF6PMg" value="48.4484" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirPMg" value="0.856" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitAMP" value="0.000255" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitAMPFDP" value="690" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitF6P" value="0.304" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitF6PMg" value="315" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitFDP" value="0.043101" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitFDPMg" value="0.00642" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitFDPMgMg" value="100" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitP" value="0.642" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitPF6P" value="0.00689" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitPF6PMg" value="16.5" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitPMg" value="539" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KmrFDP" value="0.06361410000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KmrMg" value="0.039039" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KmtFDP" value="1e-005" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KmtMg" value="55.055" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=L0" value="0.0008150000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=Vmax" value="0.215583" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=n" value="4" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KdFDPMg" value="5.81" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdFDPMg],Reference=InitialValue&gt;
              </InitialExpression>
            </ModelParameter>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KdrOAA" value="4.35404" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KdrPEP" value="655.591" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KdtOAA" value="17.9127" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KdtPEP" value="0.0121991" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KefrACCOA" value="0.140029" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KefrASP" value="0.389183" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KefrCIT" value="34.4277" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KefrCYS" value="0.000448843" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KefrFDP" value="9.991260000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KefrFDPACCOA" value="0.0156251" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KefrFUM" value="2.7475" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KefrMAL" value="0.23002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KefrSUC" value="22.9834" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KeftACCOA" value="1.27598" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KeftASP" value="27.4729" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KeftCIT" value="0.521945" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KeftCYS" value="0.977374" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KeftFDP" value="13.2064" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KeftFDPACCOA" value="47.7563" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KeftFUM" value="9.761189999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KeftMAL" value="0.737283" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KeftSUC" value="107.18" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=Keq" value="149.705" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KmrHCO3" value="0.00219811" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KmrOAA" value="13.0469" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KmrP" value="0.663356" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KmrPEP" value="3.20089" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KmtHCO3" value="0.00220031" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KmtOAA" value="6.80995" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KmtP" value="0.285131" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=KmtPEP" value="5.12497" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=L0" value="6.37209e-006" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=Vmax" value="21.439" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPC],ParameterGroup=Parameters,Parameter=n" value="4.00319" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PCK]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PCK],ParameterGroup=Parameters,Parameter=Keq" value="1.88166" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PCK],ParameterGroup=Parameters,Parameter=KmADP" value="0.04997" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PCK],ParameterGroup=Parameters,Parameter=KmATP" value="0.0600418" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PCK],ParameterGroup=Parameters,Parameter=KmHCO3" value="2.6319" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PCK],ParameterGroup=Parameters,Parameter=KmOAA" value="0.66966" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PCK],ParameterGroup=Parameters,Parameter=KmPEP" value="0.0699833" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PCK],ParameterGroup=Parameters,Parameter=Vmax" value="8.087770000000001" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdAMP" value="1480" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdATPMgPPS" value="0.0549" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdMg" value="36.9" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdP" value="346" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdPEP" value="95.7" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdPYR" value="2740" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KefADP" value="0.0283" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KefAKG" value="0.274" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KefATP" value="0.0006279999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KefOAA" value="0.796" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=Keq" value="200000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KmAMP" value="0.000384" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KmATPMg" value="0.0549" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KmP" value="84.40000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KmPEP" value="20.7" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KmPYR" value="0.229" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=Vmax" value="0.0163772" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=W" value="10" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=alpha" value="38900" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdADPMg" value="1.27771" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdADPMg],Reference=InitialValue&gt;
              </InitialExpression>
            </ModelParameter>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdATPMg" value="0.0847634" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdATPMg],Reference=InitialValue&gt;
              </InitialExpression>
            </ModelParameter>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KefrACCOA" value="1.83336" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KefrASP" value="0.362142" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KefrATP" value="88.9752" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KefrCOA" value="0.268" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KeftACCOA" value="0.197" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KeftASP" value="0.583" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KeftATP" value="0.26" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KeftCOA" value="0.268" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KirNAD" value="0.6364570000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KitNAD" value="0.990398" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KmrMAL" value="0.212913" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KmrMg" value="0.191871" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KmrMn" value="0.272568" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KmrNAD" value="1.36636" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KmtMAL" value="0.093" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KmtMg" value="2.37681" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KmtMn" value="0.410198" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=KmtNAD" value="0.108" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=L0" value="19.9" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=Vmax" value="6.64269" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MAE],ParameterGroup=Parameters,Parameter=n" value="4" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PDH]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PDH],ParameterGroup=Parameters,Parameter=Keq" value="3138.16" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PDH],ParameterGroup=Parameters,Parameter=KmACCOA" value="10.174" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PDH],ParameterGroup=Parameters,Parameter=KmCOA" value="0.005004610000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PDH],ParameterGroup=Parameters,Parameter=KmHCO3" value="0.00545112" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PDH],ParameterGroup=Parameters,Parameter=KmNAD" value="0.009990000000000002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PDH],ParameterGroup=Parameters,Parameter=KmNADH" value="6.63512" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PDH],ParameterGroup=Parameters,Parameter=KmPYR" value="2" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PDH],ParameterGroup=Parameters,Parameter=Vmax" value="961.706" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=KdACCOA0" value="0.699309" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=KdcsCIT" value="7.38128" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=KdcsCOA" value="0.00174282" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=KdcsOAA" value="0.154601" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=Keq" value="83067.10000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=Ki1AKG" value="0.0150083" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=Ki1NADH" value="0.000330313" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=Ki2AKG" value="0.256252" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=Ki2NADH" value="0.0504" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=KiATP" value="0.579939" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=KmACCOA0" value="0.119974" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=KmOAA0" value="0.00123458" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=KmcsCIT" value="1.15715" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=KmcsCOA" value="9.6344e-005" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GLT],ParameterGroup=Parameters,Parameter=Vmax" value="57.0584" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACN_1]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACN_1],ParameterGroup=Parameters,Parameter=Keq" value="0.384503" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACN_1],ParameterGroup=Parameters,Parameter=KmACO" value="0.02001" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KmACO_ACN],Reference=InitialValue&gt;
              </InitialExpression>
            </ModelParameter>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACN_1],ParameterGroup=Parameters,Parameter=KmCIT" value="0.06288820000000001" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KmCIT_ACN],Reference=InitialValue&gt;
              </InitialExpression>
            </ModelParameter>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACN_1],ParameterGroup=Parameters,Parameter=KmICIT" value="9.313520000000001" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KmICIT_ACN],Reference=InitialValue&gt;
              </InitialExpression>
            </ModelParameter>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACN_1],ParameterGroup=Parameters,Parameter=Vmax" value="9.724130000000001" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACN_2]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACN_2],ParameterGroup=Parameters,Parameter=Keq" value="3.49183" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACN_2],ParameterGroup=Parameters,Parameter=KmACO" value="0.02001" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KmACO_ACN],Reference=InitialValue&gt;
              </InitialExpression>
            </ModelParameter>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACN_2],ParameterGroup=Parameters,Parameter=KmCIT" value="0.06288820000000001" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KmCIT_ACN],Reference=InitialValue&gt;
              </InitialExpression>
            </ModelParameter>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACN_2],ParameterGroup=Parameters,Parameter=KmICIT" value="9.313520000000001" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KmICIT_ACN],Reference=InitialValue&gt;
              </InitialExpression>
            </ModelParameter>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACN_2],ParameterGroup=Parameters,Parameter=Vmax" value="9.86571" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ICD]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ICD],ParameterGroup=Parameters,Parameter=Keq" value="28.1631" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ICD],ParameterGroup=Parameters,Parameter=KmAKG" value="0.038038" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ICD],ParameterGroup=Parameters,Parameter=KmICIT" value="0.010989" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ICD],ParameterGroup=Parameters,Parameter=KmNADP" value="0.005994" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ICD],ParameterGroup=Parameters,Parameter=KmNADPH" value="0.0006833330000000002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ICD],ParameterGroup=Parameters,Parameter=kcat" value="2461.97" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[LPD]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[LPD],ParameterGroup=Parameters,Parameter=KdAKG" value="14.9386" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[LPD],ParameterGroup=Parameters,Parameter=KmAKG" value="0.020002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[LPD],ParameterGroup=Parameters,Parameter=KmCOA" value="0.0760076" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[LPD],ParameterGroup=Parameters,Parameter=KmNAD" value="0.09800979999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[LPD],ParameterGroup=Parameters,Parameter=Vmax" value="0.0684413" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[LPD],ParameterGroup=Parameters,Parameter=alpha" value="16.4304" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SK]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SK],ParameterGroup=Parameters,Parameter=Keq" value="1.15994" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SK],ParameterGroup=Parameters,Parameter=KmADP" value="0.00868475" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SK],ParameterGroup=Parameters,Parameter=KmATP" value="0.102321" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SK],ParameterGroup=Parameters,Parameter=KmCOA" value="0.255019" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SK],ParameterGroup=Parameters,Parameter=KmP" value="0.914709" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SK],ParameterGroup=Parameters,Parameter=KmSUC" value="0.800744" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SK],ParameterGroup=Parameters,Parameter=KmSUCCOA" value="0.008500000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SK],ParameterGroup=Parameters,Parameter=Vmax" value="76.8163" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SDH]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SDH],ParameterGroup=Parameters,Parameter=KefFUM" value="0.067048" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SDH],ParameterGroup=Parameters,Parameter=KefSUC" value="0.0322425" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SDH],ParameterGroup=Parameters,Parameter=Keq" value="2250" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SDH],ParameterGroup=Parameters,Parameter=KmFUM" value="1.36019" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SDH],ParameterGroup=Parameters,Parameter=KmQ" value="0.00160718" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SDH],ParameterGroup=Parameters,Parameter=KmQH2" value="0.00611663" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SDH],ParameterGroup=Parameters,Parameter=KmSUC" value="0.805727" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SDH],ParameterGroup=Parameters,Parameter=Vmax" value="1.56184" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FUMA]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FUMA],ParameterGroup=Parameters,Parameter=Keq" value="4.99966" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FUMA],ParameterGroup=Parameters,Parameter=KmFUM" value="0.6" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FUMA],ParameterGroup=Parameters,Parameter=KmMAL" value="0.7" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FUMA],ParameterGroup=Parameters,Parameter=Vmax" value="53.3414" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MQO]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MQO],ParameterGroup=Parameters,Parameter=Keq" value="9" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MQO],ParameterGroup=Parameters,Parameter=KmMAL" value="0.435" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MQO],ParameterGroup=Parameters,Parameter=KmOAA" value="75.8036" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MQO],ParameterGroup=Parameters,Parameter=KmQ" value="0.0414" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MQO],ParameterGroup=Parameters,Parameter=KmQH2" value="8.77942" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MQO],ParameterGroup=Parameters,Parameter=Vmax" value="4.62283" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MDH]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MDH],ParameterGroup=Parameters,Parameter=Keq" value="100000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MDH],ParameterGroup=Parameters,Parameter=KiNAD" value="0.0233122" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MDH],ParameterGroup=Parameters,Parameter=KiNADH" value="0.000196981" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MDH],ParameterGroup=Parameters,Parameter=KiOAA" value="2.46446" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MDH],ParameterGroup=Parameters,Parameter=KmMAL" value="0.86" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MDH],ParameterGroup=Parameters,Parameter=KmNAD" value="0.64" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MDH],ParameterGroup=Parameters,Parameter=KmNADH" value="0.003" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MDH],ParameterGroup=Parameters,Parameter=KmOAA" value="0.001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[MDH],ParameterGroup=Parameters,Parameter=Vmax" value="6.11492" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEA]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEA],ParameterGroup=Parameters,Parameter=KdICITsuc" value="0.004890739999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEA],ParameterGroup=Parameters,Parameter=KdPEP" value="1.05105" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEA],ParameterGroup=Parameters,Parameter=KdPEPglx" value="0.0358555" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEA],ParameterGroup=Parameters,Parameter=KdPEPicit" value="0.164263" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEA],ParameterGroup=Parameters,Parameter=KdPGA3" value="0.8" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEA],ParameterGroup=Parameters,Parameter=KdSUC" value="0.53053" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEA],ParameterGroup=Parameters,Parameter=Keq" value="8.8088" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEA],ParameterGroup=Parameters,Parameter=KmGLX" value="0.13013" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEA],ParameterGroup=Parameters,Parameter=KmICIT" value="0.06306299999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEA],ParameterGroup=Parameters,Parameter=KmSUC" value="0.58941" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEA],ParameterGroup=Parameters,Parameter=Vmax" value="1.52595" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEB]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEB],ParameterGroup=Parameters,Parameter=Keq" value="230000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEB],ParameterGroup=Parameters,Parameter=KmACCOA" value="0.008999999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEB],ParameterGroup=Parameters,Parameter=KmCOA" value="10.5652" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEB],ParameterGroup=Parameters,Parameter=KmGLX" value="0.021" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEB],ParameterGroup=Parameters,Parameter=KmMAL" value="15.0572" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEB],ParameterGroup=Parameters,Parameter=Vmax" value="0.352769" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEK_1]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEK_1],ParameterGroup=Parameters,Parameter=Keq" value="888" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEK_1],ParameterGroup=Parameters,Parameter=k" value="1.25457" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEK_2]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEK_2],ParameterGroup=Parameters,Parameter=Keq" value="400" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACEK_2],ParameterGroup=Parameters,Parameter=k" value="0.0332" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[EDD]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[EDD],ParameterGroup=Parameters,Parameter=Keq" value="1000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[EDD],ParameterGroup=Parameters,Parameter=KmKDPG" value="0.318316" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[EDD],ParameterGroup=Parameters,Parameter=KmPGN" value="0.6" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[EDD],ParameterGroup=Parameters,Parameter=Vmax" value="0.111359" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[EDA]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[EDA],ParameterGroup=Parameters,Parameter=Keq" value="0.5" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[EDA],ParameterGroup=Parameters,Parameter=KmGAP" value="86.678" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[EDA],ParameterGroup=Parameters,Parameter=KmKDPG" value="0.06" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[EDA],ParameterGroup=Parameters,Parameter=KmPYR" value="10" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[EDA],ParameterGroup=Parameters,Parameter=Vmax" value="0.0775241" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[NDHI]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[NDHI],ParameterGroup=Parameters,Parameter=Vmax" value="23.0735" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[NDHI],ParameterGroup=Parameters,Parameter=Keq" value="27.6193" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KeqNDH],Reference=InitialValue&gt;
              </InitialExpression>
            </ModelParameter>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PNT]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PNT],ParameterGroup=Parameters,Parameter=Keq" value="0.181891" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PNT],ParameterGroup=Parameters,Parameter=k" value="2.49441" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ADK]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ADK],ParameterGroup=Parameters,Parameter=Keq" value="0.962758" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ADK],ParameterGroup=Parameters,Parameter=k" value="0.242256" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ATP_SYN]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ATP_SYN],ParameterGroup=Parameters,Parameter=Vmax" value="108.733" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ATP_SYN],ParameterGroup=Parameters,Parameter=Keq" value="49.8315" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[CYA]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[CYA],ParameterGroup=Parameters,Parameter=Keq" value="2591.19" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[CYA],ParameterGroup=Parameters,Parameter=k" value="0.00414418" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[CYA],ParameterGroup=Parameters,Parameter=KaeiiaP" value="0.180981" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[DOS]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[DOS],ParameterGroup=Parameters,Parameter=k" value="0.008279999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[DOS],ParameterGroup=Parameters,Parameter=Keq" value="674.242" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACK]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACK],ParameterGroup=Parameters,Parameter=Keq" value="174" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACK],ParameterGroup=Parameters,Parameter=KmACE" value="7" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACK],ParameterGroup=Parameters,Parameter=KmACP" value="0.16" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACK],ParameterGroup=Parameters,Parameter=KmADP" value="0.5" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACK],ParameterGroup=Parameters,Parameter=KmATP" value="0.07000000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACK],ParameterGroup=Parameters,Parameter=Vmax" value="7.23" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACS]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACS],ParameterGroup=Parameters,Parameter=KmACE" value="0.07000000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACS],ParameterGroup=Parameters,Parameter=KmATP" value="0.1" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACS],ParameterGroup=Parameters,Parameter=KmCOA" value="0.01" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ACS],ParameterGroup=Parameters,Parameter=Vmax" value="7.3" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTA]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTA],ParameterGroup=Parameters,Parameter=Keq" value="0.0053952" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTA],ParameterGroup=Parameters,Parameter=KiACCOA" value="0.2" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTA],ParameterGroup=Parameters,Parameter=KiACP" value="0.2" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTA],ParameterGroup=Parameters,Parameter=KiCOA" value="0.029" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTA],ParameterGroup=Parameters,Parameter=KiP" value="13.54" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTA],ParameterGroup=Parameters,Parameter=KmACP" value="0.7" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTA],ParameterGroup=Parameters,Parameter=KmP" value="6.1" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTA],ParameterGroup=Parameters,Parameter=Vmax" value="2.7" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_0]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_0],ParameterGroup=Parameters,Parameter=KmPEP" value="0.6" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_0],ParameterGroup=Parameters,Parameter=KmPYR" value="1" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_0],ParameterGroup=Parameters,Parameter=kF" value="12000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_0],ParameterGroup=Parameters,Parameter=kR" value="8000" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_1]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_1],ParameterGroup=Parameters,Parameter=k1" value="200000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_1],ParameterGroup=Parameters,Parameter=k2" value="8000" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_2]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_2],ParameterGroup=Parameters,Parameter=k1" value="61000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_2],ParameterGroup=Parameters,Parameter=k2" value="47000" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_3]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_3],ParameterGroup=Parameters,Parameter=k1" value="11000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_3],ParameterGroup=Parameters,Parameter=k2" value="4000" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_4]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_4],ParameterGroup=Parameters,Parameter=KmG6P" value="2125.91" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_4],ParameterGroup=Parameters,Parameter=KmGLC" value="0.02" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_4],ParameterGroup=Parameters,Parameter=kF" value="4000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_4],ParameterGroup=Parameters,Parameter=kR" value="1e-005" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[_GLC_FEED]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[_GLC_FEED],ParameterGroup=Parameters,Parameter=FEED" value="0.23" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[FEED],Reference=InitialValue&gt;
              </InitialExpression>
            </ModelParameter>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[CYTBO]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[CYTBO],ParameterGroup=Parameters,Parameter=Vmax" value="8.54045" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[CYTBO],ParameterGroup=Parameters,Parameter=Keq" value="12.067" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SQR]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SQR],ParameterGroup=Parameters,Parameter=Keq" value="0.94033" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[SQR],ParameterGroup=Parameters,Parameter=Vmax" value="3.41617" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[NDHII]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[NDHII],ParameterGroup=Parameters,Parameter=Vmax" value="30.8306" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[NDHII],ParameterGroup=Parameters,Parameter=Keq" value="27.6193" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KeqNDH],Reference=InitialValue&gt;
              </InitialExpression>
            </ModelParameter>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=Vmax" value="9.74137" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmG6P" value="1.20911" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmE4P" value="1.63298" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmPGA3" value="0.0764619" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmOAA" value="0.0248351" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmAKG" value="5.11989" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmPYR" value="0.00463904" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmR5P" value="0.021234" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmPEP" value="0.457734" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmGAP" value="0.024854" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmF6P" value="0.366423" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmNADPH" value="3.59774" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmACCOA" value="0.0494404" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmNAD" value="2.82239" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],ParameterGroup=Parameters,Parameter=KmATP" value="0.0468266" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ATP_MAINTENANCE]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ATP_MAINTENANCE],ParameterGroup=Parameters,Parameter=Vmax" value="1.30166" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ATP_MAINTENANCE],ParameterGroup=Parameters,Parameter=Keq" value="3.63369" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_GLC]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_GLC],ParameterGroup=Parameters,Parameter=Vmax" value="100" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_GLC],ParameterGroup=Parameters,Parameter=Km" value="10" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PIT]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PIT],ParameterGroup=Parameters,Parameter=Vmax" value="7.146" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PIT],ParameterGroup=Parameters,Parameter=KmP" value="12.18" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PIT],ParameterGroup=Parameters,Parameter=KmPp" value="0.025" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PIT],ParameterGroup=Parameters,Parameter=Kr" value="0.1" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_P]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_P],ParameterGroup=Parameters,Parameter=P" value="1.505449345660231e+021" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_periplasm],Vector=Metabolites[Pp],Reference=InitialParticleNumber&gt;
              </InitialExpression>
            </ModelParameter>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_P],ParameterGroup=Parameters,Parameter=Vmax" value="100" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_P],ParameterGroup=Parameters,Parameter=Km" value="10" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_ACE1]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_ACE1],ParameterGroup=Parameters,Parameter=P" value="2.810796198627631e+016" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_periplasm],Vector=Metabolites[ACEp],Reference=InitialParticleNumber&gt;
              </InitialExpression>
            </ModelParameter>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_ACE1],ParameterGroup=Parameters,Parameter=Vmax" value="100" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_ACE1],ParameterGroup=Parameters,Parameter=Km" value="10" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[_ACE_OUT]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[_ACE_OUT],ParameterGroup=Parameters,Parameter=k1" value="5.556e-005" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_ACE2]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_ACE2],ParameterGroup=Parameters,Parameter=P" value="1.124256015786816e+019" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[extracellular],Vector=Metabolites[ACEx],Reference=InitialParticleNumber&gt;
              </InitialExpression>
            </ModelParameter>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_ACE2],ParameterGroup=Parameters,Parameter=Vmax" value="100" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_ACE2],ParameterGroup=Parameters,Parameter=Km" value="10" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GL6P_HYDROLYSIS]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GL6P_HYDROLYSIS],ParameterGroup=Parameters,Parameter=KGl6Phydrol" value="0.000167" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GL6P_HYDROLYSIS],ParameterGroup=Parameters,Parameter=KeqGl6Phydrol" value="42.8" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
        </ModelParameterGroup>
      </ModelParameterSet>
    </ListOfModelParameterSets>
    <StateTemplate>
      <StateTemplateVariable objectReference="Model_3"/>
      <StateTemplateVariable objectReference="Metabolite_113"/>
      <StateTemplateVariable objectReference="Metabolite_147"/>
      <StateTemplateVariable objectReference="Metabolite_101"/>
      <StateTemplateVariable objectReference="Metabolite_57"/>
      <StateTemplateVariable objectReference="Metabolite_37"/>
      <StateTemplateVariable objectReference="Metabolite_55"/>
      <StateTemplateVariable objectReference="Metabolite_25"/>
      <StateTemplateVariable objectReference="Metabolite_41"/>
      <StateTemplateVariable objectReference="Metabolite_17"/>
      <StateTemplateVariable objectReference="Metabolite_1"/>
      <StateTemplateVariable objectReference="Metabolite_47"/>
      <StateTemplateVariable objectReference="Metabolite_53"/>
      <StateTemplateVariable objectReference="Metabolite_35"/>
      <StateTemplateVariable objectReference="Metabolite_67"/>
      <StateTemplateVariable objectReference="Metabolite_45"/>
      <StateTemplateVariable objectReference="Metabolite_139"/>
      <StateTemplateVariable objectReference="Metabolite_97"/>
      <StateTemplateVariable objectReference="Metabolite_63"/>
      <StateTemplateVariable objectReference="Metabolite_23"/>
      <StateTemplateVariable objectReference="Metabolite_103"/>
      <StateTemplateVariable objectReference="Metabolite_31"/>
      <StateTemplateVariable objectReference="Metabolite_13"/>
      <StateTemplateVariable objectReference="Metabolite_77"/>
      <StateTemplateVariable objectReference="Metabolite_141"/>
      <StateTemplateVariable objectReference="Metabolite_145"/>
      <StateTemplateVariable objectReference="Metabolite_151"/>
      <StateTemplateVariable objectReference="Metabolite_51"/>
      <StateTemplateVariable objectReference="Metabolite_73"/>
      <StateTemplateVariable objectReference="Metabolite_89"/>
      <StateTemplateVariable objectReference="Metabolite_65"/>
      <StateTemplateVariable objectReference="Metabolite_3"/>
      <StateTemplateVariable objectReference="Metabolite_27"/>
      <StateTemplateVariable objectReference="Metabolite_69"/>
      <StateTemplateVariable objectReference="Metabolite_135"/>
      <StateTemplateVariable objectReference="Metabolite_15"/>
      <StateTemplateVariable objectReference="Metabolite_49"/>
      <StateTemplateVariable objectReference="Metabolite_33"/>
      <StateTemplateVariable objectReference="Metabolite_153"/>
      <StateTemplateVariable objectReference="Metabolite_21"/>
      <StateTemplateVariable objectReference="Metabolite_81"/>
      <StateTemplateVariable objectReference="Metabolite_71"/>
      <StateTemplateVariable objectReference="Metabolite_5"/>
      <StateTemplateVariable objectReference="Metabolite_9"/>
      <StateTemplateVariable objectReference="Metabolite_149"/>
      <StateTemplateVariable objectReference="Metabolite_29"/>
      <StateTemplateVariable objectReference="Metabolite_61"/>
      <StateTemplateVariable objectReference="Metabolite_87"/>
      <StateTemplateVariable objectReference="Metabolite_95"/>
      <StateTemplateVariable objectReference="Metabolite_11"/>
      <StateTemplateVariable objectReference="Metabolite_107"/>
      <StateTemplateVariable objectReference="Metabolite_7"/>
      <StateTemplateVariable objectReference="Metabolite_19"/>
      <StateTemplateVariable objectReference="Metabolite_105"/>
      <StateTemplateVariable objectReference="Metabolite_91"/>
      <StateTemplateVariable objectReference="Metabolite_43"/>
      <StateTemplateVariable objectReference="Metabolite_79"/>
      <StateTemplateVariable objectReference="Metabolite_83"/>
      <StateTemplateVariable objectReference="Metabolite_137"/>
      <StateTemplateVariable objectReference="Metabolite_93"/>
      <StateTemplateVariable objectReference="Metabolite_75"/>
      <StateTemplateVariable objectReference="Metabolite_85"/>
      <StateTemplateVariable objectReference="Metabolite_59"/>
      <StateTemplateVariable objectReference="Metabolite_39"/>
      <StateTemplateVariable objectReference="Metabolite_99"/>
      <StateTemplateVariable objectReference="Metabolite_117"/>
      <StateTemplateVariable objectReference="Metabolite_119"/>
      <StateTemplateVariable objectReference="Metabolite_121"/>
      <StateTemplateVariable objectReference="Metabolite_109"/>
      <StateTemplateVariable objectReference="Metabolite_111"/>
      <StateTemplateVariable objectReference="Metabolite_115"/>
      <StateTemplateVariable objectReference="Metabolite_123"/>
      <StateTemplateVariable objectReference="Metabolite_125"/>
      <StateTemplateVariable objectReference="Metabolite_127"/>
      <StateTemplateVariable objectReference="Metabolite_129"/>
      <StateTemplateVariable objectReference="Metabolite_131"/>
      <StateTemplateVariable objectReference="Metabolite_133"/>
      <StateTemplateVariable objectReference="Metabolite_143"/>
      <StateTemplateVariable objectReference="Compartment_1"/>
      <StateTemplateVariable objectReference="Compartment_3"/>
      <StateTemplateVariable objectReference="Compartment_5"/>
      <StateTemplateVariable objectReference="ModelValue_0"/>
      <StateTemplateVariable objectReference="ModelValue_1"/>
      <StateTemplateVariable objectReference="ModelValue_2"/>
      <StateTemplateVariable objectReference="ModelValue_3"/>
      <StateTemplateVariable objectReference="ModelValue_4"/>
      <StateTemplateVariable objectReference="ModelValue_5"/>
      <StateTemplateVariable objectReference="ModelValue_6"/>
      <StateTemplateVariable objectReference="ModelValue_7"/>
    </StateTemplate>
    <InitialState type="initialState">
      0 5.879987842207703e+021 8091473706854326 3.603135976315973e+020 1.963925127694926e+020 8.500516016790164e+020 1.426593487227837e+020 7.05693301383649e+019 1.010739458841187e+020 1.576394850846697e+020 9.318843183080258e+019 6.004306279627797e+020 7.925125419323508e+019 6.215772458634581e+020 1.301115620269875e+020 7.698681357565698e+019 1.124380943115288e+017 8.329093572866015e+018 2.058528180838142e+020 5.185843963253621e+020 1.121639644055652e+020 6.370367694157086e+019 2.632242451004475e+020 8.552558885904762e+018 3.814045050430201e+020 1.124256015786816e+019 1.505449345660231e+021 4.193061747813792e+020 2.011473853979885e+017 3.268480204715989e+017 8.550546593595767e+019 1.936106650655732e+019 1.964213101596648e+018 2.474364250168681e+019 1.529107897748966e+020 7.907068085616183e+019 2.278158803291449e+020 5.228284334351698e+019 2.810796198627631e+016 1.283408617631331e+020 2.841906957912192e+016 3.047312090476589e+020 8.666420504823286e+018 3.939117017287109e+019 6.072381096825504e+017 5.305754012188268e+018 6.434157914056002e+019 3.162113325697529e+018 3.445812690814797e+019 5.39165207053416e+019 5.55922809558007e+020 3.600459862224256e+020 1.697088696281518e+020 1.549016191114563e+021 3.098828960171077e+019 5.369509789454854e+019 2.894051911487251e+020 2.063314616344276e+017 4.493033892713735e+020 1.674723810868397e+018 3.842826679553765e+018 1.15150454342188e+017 4.058216660206235e+020 9.542465977025343e+019 3.382589890806846e+019 1.581911646485274e+020 1.42797608318511e+021 2.492053885875944e+019 3.011070895e+020 8.430998506e+020 6.022141789999999e+020 7.045905894300009e+020 5.118820521500001e+019 1.806642537000001e+020 1.904368444786942e+016 6.022141789999999e+020 1.2646497759e+020 6.022141790000004e+023 1 100 0.25 0.23 1.27771 0.0847634 5.81 9.313520000000001 0.06288820000000001 0.02001 27.6193 
    </InitialState>
  </Model>
  <ListOfTasks>
    <Task key="Task_14" name="Steady-State" type="steadyState" scheduled="true" updateModel="false">
      <Report reference="Report_9" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="JacobianRequested" type="bool" value="1"/>
        <Parameter name="StabilityAnalysisRequested" type="bool" value="1"/>
      </Problem>
      <Method name="Enhanced Newton" type="EnhancedNewton">
        <Parameter name="Resolution" type="unsignedFloat" value="1e-009"/>
        <Parameter name="Derivation Factor" type="unsignedFloat" value="0.001"/>
        <Parameter name="Use Newton" type="bool" value="0"/>
        <Parameter name="Use Integration" type="bool" value="1"/>
        <Parameter name="Use Back Integration" type="bool" value="1"/>
        <Parameter name="Accept Negative Concentrations" type="bool" value="0"/>
        <Parameter name="Iteration Limit" type="unsignedInteger" value="50"/>
        <Parameter name="Maximum duration for forward integration" type="unsignedFloat" value="1000000000"/>
        <Parameter name="Maximum duration for backward integration" type="unsignedFloat" value="1000000"/>
      </Method>
    </Task>
    <Task key="Task_15" name="Time-Course" type="timeCourse" scheduled="false" updateModel="false">
      <Problem>
        <Parameter name="StepNumber" type="unsignedInteger" value="10000"/>
        <Parameter name="StepSize" type="float" value="0.05"/>
        <Parameter name="Duration" type="float" value="500"/>
        <Parameter name="TimeSeriesRequested" type="bool" value="0"/>
        <Parameter name="OutputStartTime" type="float" value="0"/>
        <Parameter name="Output Event" type="bool" value="0"/>
        <Parameter name="Continue on Simultaneous Events" type="bool" value="0"/>
      </Problem>
      <Method name="Deterministic (LSODA)" type="Deterministic(LSODA)">
        <Parameter name="Integrate Reduced Model" type="bool" value="0"/>
        <Parameter name="Relative Tolerance" type="unsignedFloat" value="1e-009"/>
        <Parameter name="Absolute Tolerance" type="unsignedFloat" value="1e-012"/>
        <Parameter name="Max Internal Steps" type="unsignedInteger" value="50000"/>
      </Method>
    </Task>
    <Task key="Task_16" name="Scan" type="scan" scheduled="false" updateModel="false">
      <Problem>
        <Parameter name="Subtask" type="unsignedInteger" value="0"/>
        <ParameterGroup name="ScanItems">
        </ParameterGroup>
        <Parameter name="Output in subtask" type="bool" value="0"/>
        <Parameter name="Adjust initial conditions" type="bool" value="0"/>
      </Problem>
      <Method name="Scan Framework" type="ScanFramework">
      </Method>
    </Task>
    <Task key="Task_17" name="Elementary Flux Modes" type="fluxMode" scheduled="false" updateModel="false">
      <Report reference="Report_10" target="" append="1" confirmOverwrite="1"/>
      <Problem>
      </Problem>
      <Method name="EFM Algorithm" type="EFMAlgorithm">
      </Method>
    </Task>
    <Task key="Task_18" name="Optimization" type="optimization" scheduled="false" updateModel="false">
      <Report reference="Report_11" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="Subtask" type="cn" value="CN=Root,Vector=TaskList[Steady-State]"/>
        <ParameterText name="ObjectiveExpression" type="expression">
          
        </ParameterText>
        <Parameter name="Maximize" type="bool" value="0"/>
        <Parameter name="Randomize Start Values" type="bool" value="0"/>
        <Parameter name="Calculate Statistics" type="bool" value="1"/>
        <ParameterGroup name="OptimizationItemList">
        </ParameterGroup>
        <ParameterGroup name="OptimizationConstraintList">
        </ParameterGroup>
      </Problem>
      <Method name="Random Search" type="RandomSearch">
        <Parameter name="Number of Iterations" type="unsignedInteger" value="100000"/>
        <Parameter name="Random Number Generator" type="unsignedInteger" value="1"/>
        <Parameter name="Seed" type="unsignedInteger" value="0"/>
      </Method>
    </Task>
    <Task key="Task_19" name="Parameter Estimation" type="parameterFitting" scheduled="false" updateModel="false">
      <Report reference="Report_12" target="31.txt" append="0" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="Maximize" type="bool" value="0"/>
        <Parameter name="Randomize Start Values" type="bool" value="0"/>
        <Parameter name="Calculate Statistics" type="bool" value="0"/>
        <ParameterGroup name="OptimizationItemList">
        </ParameterGroup>
        <ParameterGroup name="OptimizationConstraintList">
        </ParameterGroup>
        <Parameter name="Steady-State" type="cn" value="CN=Root,Vector=TaskList[Steady-State]"/>
        <Parameter name="Time-Course" type="cn" value="CN=Root,Vector=TaskList[Time-Course]"/>
        <Parameter name="Create Parameter Sets" type="bool" value="0"/>
        <ParameterGroup name="Experiment Set">
        </ParameterGroup>
        <ParameterGroup name="Validation Set">
          <Parameter name="Threshold" type="unsignedInteger" value="5"/>
          <Parameter name="Weight" type="unsignedFloat" value="1"/>
        </ParameterGroup>
      </Problem>
      <Method name="Hooke &amp; Jeeves" type="HookeJeeves">
        <Parameter name="Iteration Limit" type="unsignedInteger" value="50"/>
        <Parameter name="Tolerance" type="float" value="1e-005"/>
        <Parameter name="Rho" type="float" value="0.2"/>
      </Method>
    </Task>
    <Task key="Task_20" name="Metabolic Control Analysis" type="metabolicControlAnalysis" scheduled="false" updateModel="false">
      <Report reference="Report_13" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="Steady-State" type="key" value="Task_14"/>
      </Problem>
      <Method name="MCA Method (Reder)" type="MCAMethod(Reder)">
        <Parameter name="Modulation Factor" type="unsignedFloat" value="1e-009"/>
      </Method>
    </Task>
    <Task key="Task_21" name="Lyapunov Exponents" type="lyapunovExponents" scheduled="false" updateModel="false">
      <Report reference="Report_14" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="ExponentNumber" type="unsignedInteger" value="3"/>
        <Parameter name="DivergenceRequested" type="bool" value="1"/>
        <Parameter name="TransientTime" type="float" value="0"/>
      </Problem>
      <Method name="Wolf Method" type="WolfMethod">
        <Parameter name="Orthonormalization Interval" type="unsignedFloat" value="1"/>
        <Parameter name="Overall time" type="unsignedFloat" value="1000"/>
        <Parameter name="Relative Tolerance" type="unsignedFloat" value="1e-006"/>
        <Parameter name="Absolute Tolerance" type="unsignedFloat" value="1e-012"/>
        <Parameter name="Max Internal Steps" type="unsignedInteger" value="10000"/>
      </Method>
    </Task>
    <Task key="Task_22" name="Time Scale Separation Analysis" type="timeScaleSeparationAnalysis" scheduled="false" updateModel="false">
      <Report reference="Report_15" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="StepNumber" type="unsignedInteger" value="100"/>
        <Parameter name="StepSize" type="float" value="0.01"/>
        <Parameter name="Duration" type="float" value="1"/>
        <Parameter name="TimeSeriesRequested" type="bool" value="1"/>
        <Parameter name="OutputStartTime" type="float" value="0"/>
      </Problem>
      <Method name="ILDM (LSODA,Deuflhard)" type="TimeScaleSeparation(ILDM,Deuflhard)">
        <Parameter name="Deuflhard Tolerance" type="unsignedFloat" value="1e-006"/>
      </Method>
    </Task>
    <Task key="Task_23" name="Sensitivities" type="sensitivities" scheduled="false" updateModel="false">
      <Report reference="Report_16" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="SubtaskType" type="unsignedInteger" value="1"/>
        <ParameterGroup name="TargetFunctions">
          <Parameter name="SingleObject" type="cn" value=""/>
          <Parameter name="ObjectListType" type="unsignedInteger" value="7"/>
        </ParameterGroup>
        <ParameterGroup name="ListOfVariables">
          <ParameterGroup name="Variables">
            <Parameter name="SingleObject" type="cn" value=""/>
            <Parameter name="ObjectListType" type="unsignedInteger" value="41"/>
          </ParameterGroup>
        </ParameterGroup>
      </Problem>
      <Method name="Sensitivities Method" type="SensitivitiesMethod">
        <Parameter name="Delta factor" type="unsignedFloat" value="0.001"/>
        <Parameter name="Delta minimum" type="unsignedFloat" value="1e-012"/>
      </Method>
    </Task>
    <Task key="Task_24" name="Moieties" type="moieties" scheduled="false" updateModel="false">
      <Problem>
      </Problem>
      <Method name="Householder Reduction" type="Householder">
      </Method>
    </Task>
    <Task key="Task_25" name="Cross Section" type="crosssection" scheduled="false" updateModel="false">
      <Problem>
        <Parameter name="StepNumber" type="unsignedInteger" value="100"/>
        <Parameter name="StepSize" type="float" value="0.01"/>
        <Parameter name="Duration" type="float" value="1"/>
        <Parameter name="TimeSeriesRequested" type="bool" value="1"/>
        <Parameter name="OutputStartTime" type="float" value="0"/>
        <Parameter name="Output Event" type="bool" value="0"/>
        <Parameter name="Continue on Simultaneous Events" type="bool" value="0"/>
        <Parameter name="LimitCrossings" type="bool" value="0"/>
        <Parameter name="NumCrossingsLimit" type="unsignedInteger" value="0"/>
        <Parameter name="LimitOutTime" type="bool" value="0"/>
        <Parameter name="LimitOutCrossings" type="bool" value="0"/>
        <Parameter name="PositiveDirection" type="bool" value="1"/>
        <Parameter name="NumOutCrossingsLimit" type="unsignedInteger" value="0"/>
        <Parameter name="LimitUntilConvergence" type="bool" value="0"/>
        <Parameter name="ConvergenceTolerance" type="float" value="0"/>
        <Parameter name="Threshold" type="float" value="0"/>
        <Parameter name="DelayOutputUntilConvergence" type="bool" value="0"/>
        <Parameter name="OutputConvergenceTolerance" type="float" value="0"/>
        <ParameterText name="TriggerExpression" type="expression">
          
        </ParameterText>
        <Parameter name="SingleVariable" type="cn" value="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GROWTH],Reference=Flux"/>
      </Problem>
      <Method name="Deterministic (LSODA)" type="Deterministic(LSODA)">
        <Parameter name="Integrate Reduced Model" type="bool" value="0"/>
        <Parameter name="Relative Tolerance" type="unsignedFloat" value="1e-006"/>
        <Parameter name="Absolute Tolerance" type="unsignedFloat" value="1e-012"/>
        <Parameter name="Max Internal Steps" type="unsignedInteger" value="10000"/>
      </Method>
    </Task>
    <Task key="Task_26" name="Linear Noise Approximation" type="linearNoiseApproximation" scheduled="false" updateModel="false">
      <Report reference="Report_17" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="Steady-State" type="key" value="Task_14"/>
      </Problem>
      <Method name="Linear Noise Approximation" type="LinearNoiseApproximation">
      </Method>
    </Task>
  </ListOfTasks>
  <ListOfReports>
    <Report key="Report_9" name="Steady-State" taskType="steadyState" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Footer>
        <Object cn="CN=Root,Vector=TaskList[Steady-State]"/>
      </Footer>
    </Report>
    <Report key="Report_10" name="Elementary Flux Modes" taskType="fluxMode" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Footer>
        <Object cn="CN=Root,Vector=TaskList[Elementary Flux Modes],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_11" name="Optimization" taskType="optimization" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Optimization],Object=Description"/>
        <Object cn="String=\[Function Evaluations\]"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="String=\[Best Value\]"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="String=\[Best Parameters\]"/>
      </Header>
      <Body>
        <Object cn="CN=Root,Vector=TaskList[Optimization],Problem=Optimization,Reference=Function Evaluations"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="CN=Root,Vector=TaskList[Optimization],Problem=Optimization,Reference=Best Value"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="CN=Root,Vector=TaskList[Optimization],Problem=Optimization,Reference=Best Parameters"/>
      </Body>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Optimization],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_12" name="Parameter Estimation" taskType="parameterFitting" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Parameter Estimation],Object=Description"/>
        <Object cn="String=\[Function Evaluations\]"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="String=\[Best Value\]"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="String=\[Best Parameters\]"/>
      </Header>
      <Body>
        <Object cn="CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,Reference=Function Evaluations"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,Reference=Best Value"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,Reference=Best Parameters"/>
      </Body>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Parameter Estimation],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_13" name="Metabolic Control Analysis" taskType="metabolicControlAnalysis" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Metabolic Control Analysis],Object=Description"/>
      </Header>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Metabolic Control Analysis],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_14" name="Lyapunov Exponents" taskType="lyapunovExponents" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Lyapunov Exponents],Object=Description"/>
      </Header>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Lyapunov Exponents],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_15" name="Time Scale Separation Analysis" taskType="timeScaleSeparationAnalysis" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Time Scale Separation Analysis],Object=Description"/>
      </Header>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Time Scale Separation Analysis],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_16" name="Sensitivities" taskType="sensitivities" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Sensitivities],Object=Description"/>
      </Header>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Sensitivities],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_17" name="Linear Noise Approximation" taskType="linearNoiseApproximation" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Linear Noise Approximation],Object=Description"/>
      </Header>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Linear Noise Approximation],Object=Result"/>
      </Footer>
    </Report>
  </ListOfReports>
  <GUI>
  </GUI>
  <SBMLReference file="ThEcoli_BIOMODELS.xml">
    <SBMLMap SBMLid="ACCOA" COPASIkey="Metabolite_1"/>
    <SBMLMap SBMLid="ACEA" COPASIkey="Reaction_35"/>
    <SBMLMap SBMLid="ACEB" COPASIkey="Reaction_36"/>
    <SBMLMap SBMLid="ACEK_1" COPASIkey="Reaction_37"/>
    <SBMLMap SBMLid="ACEK_2" COPASIkey="Reaction_38"/>
    <SBMLMap SBMLid="ACEp" COPASIkey="Metabolite_153"/>
    <SBMLMap SBMLid="ACEx" COPASIkey="Metabolite_139"/>
    <SBMLMap SBMLid="ACEx_0" COPASIkey="Metabolite_145"/>
    <SBMLMap SBMLid="ACK" COPASIkey="Reaction_47"/>
    <SBMLMap SBMLid="ACN_1" COPASIkey="Reaction_26"/>
    <SBMLMap SBMLid="ACN_2" COPASIkey="Reaction_27"/>
    <SBMLMap SBMLid="ACO" COPASIkey="Metabolite_3"/>
    <SBMLMap SBMLid="ACP" COPASIkey="Metabolite_5"/>
    <SBMLMap SBMLid="ACS" COPASIkey="Reaction_48"/>
    <SBMLMap SBMLid="ADK" COPASIkey="Reaction_43"/>
    <SBMLMap SBMLid="ADP" COPASIkey="Metabolite_101"/>
    <SBMLMap SBMLid="AKG" COPASIkey="Metabolite_7"/>
    <SBMLMap SBMLid="AMP" COPASIkey="Metabolite_103"/>
    <SBMLMap SBMLid="ASP" COPASIkey="Metabolite_123"/>
    <SBMLMap SBMLid="ATP" COPASIkey="Metabolite_105"/>
    <SBMLMap SBMLid="ATP_MAINTENANCE" COPASIkey="Reaction_60"/>
    <SBMLMap SBMLid="ATP_syn" COPASIkey="Reaction_44"/>
    <SBMLMap SBMLid="BPG" COPASIkey="Metabolite_9"/>
    <SBMLMap SBMLid="CAMP" COPASIkey="Metabolite_107"/>
    <SBMLMap SBMLid="CIT" COPASIkey="Metabolite_11"/>
    <SBMLMap SBMLid="COA" COPASIkey="Metabolite_109"/>
    <SBMLMap SBMLid="CYA" COPASIkey="Reaction_45"/>
    <SBMLMap SBMLid="CYS" COPASIkey="Metabolite_125"/>
    <SBMLMap SBMLid="CYTBO" COPASIkey="Reaction_56"/>
    <SBMLMap SBMLid="DAP" COPASIkey="Metabolite_13"/>
    <SBMLMap SBMLid="DOS" COPASIkey="Reaction_46"/>
    <SBMLMap SBMLid="E4P" COPASIkey="Metabolite_15"/>
    <SBMLMap SBMLid="EDA" COPASIkey="Reaction_40"/>
    <SBMLMap SBMLid="EDD" COPASIkey="Reaction_39"/>
    <SBMLMap SBMLid="ENO" COPASIkey="Reaction_7"/>
    <SBMLMap SBMLid="F6P" COPASIkey="Metabolite_17"/>
    <SBMLMap SBMLid="F6P_E4P_TKT" COPASIkey="Reaction_15"/>
    <SBMLMap SBMLid="F6P_GAP_TAL" COPASIkey="Reaction_17"/>
    <SBMLMap SBMLid="FAD" COPASIkey="Metabolite_135"/>
    <SBMLMap SBMLid="FADH2" COPASIkey="Metabolite_137"/>
    <SBMLMap SBMLid="FBA" COPASIkey="Reaction_2"/>
    <SBMLMap SBMLid="FBP" COPASIkey="Reaction_19"/>
    <SBMLMap SBMLid="FDP" COPASIkey="Metabolite_19"/>
    <SBMLMap SBMLid="FEED" COPASIkey="ModelValue_0"/>
    <SBMLMap SBMLid="FUM" COPASIkey="Metabolite_21"/>
    <SBMLMap SBMLid="FUMA" COPASIkey="Reaction_32"/>
    <SBMLMap SBMLid="Function_for_ACEA" COPASIkey="Function_75"/>
    <SBMLMap SBMLid="Function_for_ACEB" COPASIkey="Function_76"/>
    <SBMLMap SBMLid="Function_for_ACEK_1" COPASIkey="Function_77"/>
    <SBMLMap SBMLid="Function_for_ACEK_2" COPASIkey="Function_78"/>
    <SBMLMap SBMLid="Function_for_ACK" COPASIkey="Function_85"/>
    <SBMLMap SBMLid="Function_for_ACN_1" COPASIkey="Function_66"/>
    <SBMLMap SBMLid="Function_for_ACN_2" COPASIkey="Function_67"/>
    <SBMLMap SBMLid="Function_for_ACS" COPASIkey="Function_86"/>
    <SBMLMap SBMLid="Function_for_ADK" COPASIkey="Function_82"/>
    <SBMLMap SBMLid="Function_for_ATP_MAINTENANCE" COPASIkey="Function_91"/>
    <SBMLMap SBMLid="Function_for_ATP_SYN" COPASIkey="Function_92"/>
    <SBMLMap SBMLid="Function_for_CYA" COPASIkey="Function_83"/>
    <SBMLMap SBMLid="Function_for_CYTBO" COPASIkey="Function_94"/>
    <SBMLMap SBMLid="Function_for_DOS" COPASIkey="Function_84"/>
    <SBMLMap SBMLid="Function_for_EDA" COPASIkey="Function_80"/>
    <SBMLMap SBMLid="Function_for_EDD" COPASIkey="Function_79"/>
    <SBMLMap SBMLid="Function_for_ENO" COPASIkey="Function_47"/>
    <SBMLMap SBMLid="Function_for_F6P_E4P_TKT" COPASIkey="Function_55"/>
    <SBMLMap SBMLid="Function_for_F6P_GAP_TAL" COPASIkey="Function_57"/>
    <SBMLMap SBMLid="Function_for_FBA" COPASIkey="Function_42"/>
    <SBMLMap SBMLid="Function_for_FBP" COPASIkey="Function_59"/>
    <SBMLMap SBMLid="Function_for_FUMA" COPASIkey="Function_72"/>
    <SBMLMap SBMLid="Function_for_GDH" COPASIkey="Function_44"/>
    <SBMLMap SBMLid="Function_for_GL6P_HYDROLYSIS" COPASIkey="Function_99"/>
    <SBMLMap SBMLid="Function_for_GLC_FEED" COPASIkey="Function_90"/>
    <SBMLMap SBMLid="Function_for_GLT" COPASIkey="Function_65"/>
    <SBMLMap SBMLid="Function_for_GND" COPASIkey="Function_51"/>
    <SBMLMap SBMLid="Function_for_GPM" COPASIkey="Function_46"/>
    <SBMLMap SBMLid="Function_for_GROWTH" COPASIkey="Function_97"/>
    <SBMLMap SBMLid="Function_for_ICD" COPASIkey="Function_68"/>
    <SBMLMap SBMLid="Function_for_LPD" COPASIkey="Function_69"/>
    <SBMLMap SBMLid="Function_for_MAD" COPASIkey="Function_63"/>
    <SBMLMap SBMLid="Function_for_MDH" COPASIkey="Function_74"/>
    <SBMLMap SBMLid="Function_for_MQO" COPASIkey="Function_73"/>
    <SBMLMap SBMLid="Function_for_NDH1" COPASIkey="Function_93"/>
    <SBMLMap SBMLid="Function_for_NDH2" COPASIkey="Function_96"/>
    <SBMLMap SBMLid="Function_for_PCK" COPASIkey="Function_61"/>
    <SBMLMap SBMLid="Function_for_PDH" COPASIkey="Function_64"/>
    <SBMLMap SBMLid="Function_for_PFK" COPASIkey="Function_41"/>
    <SBMLMap SBMLid="Function_for_PGI" COPASIkey="Function_40"/>
    <SBMLMap SBMLid="Function_for_PGK" COPASIkey="Function_45"/>
    <SBMLMap SBMLid="Function_for_PGL" COPASIkey="Function_50"/>
    <SBMLMap SBMLid="Function_for_PIT" COPASIkey="Function_98"/>
    <SBMLMap SBMLid="Function_for_PNT" COPASIkey="Function_81"/>
    <SBMLMap SBMLid="Function_for_PPC" COPASIkey="Function_60"/>
    <SBMLMap SBMLid="Function_for_PPS" COPASIkey="Function_62"/>
    <SBMLMap SBMLid="Function_for_PTA" COPASIkey="Function_87"/>
    <SBMLMap SBMLid="Function_for_PTS_0" COPASIkey="Function_88"/>
    <SBMLMap SBMLid="Function_for_PTS_4" COPASIkey="Function_89"/>
    <SBMLMap SBMLid="Function_for_PYK" COPASIkey="Function_48"/>
    <SBMLMap SBMLid="Function_for_RPE" COPASIkey="Function_52"/>
    <SBMLMap SBMLid="Function_for_RPI" COPASIkey="Function_53"/>
    <SBMLMap SBMLid="Function_for_S7P_E4P_TAL" COPASIkey="Function_58"/>
    <SBMLMap SBMLid="Function_for_S7P_R5P_TKT" COPASIkey="Function_56"/>
    <SBMLMap SBMLid="Function_for_SDH" COPASIkey="Function_71"/>
    <SBMLMap SBMLid="Function_for_SK" COPASIkey="Function_70"/>
    <SBMLMap SBMLid="Function_for_SQR" COPASIkey="Function_95"/>
    <SBMLMap SBMLid="Function_for_TPI" COPASIkey="Function_43"/>
    <SBMLMap SBMLid="Function_for_X5P_GAP_TKT" COPASIkey="Function_54"/>
    <SBMLMap SBMLid="Function_for_ZWF" COPASIkey="Function_49"/>
    <SBMLMap SBMLid="G6P" COPASIkey="Metabolite_23"/>
    <SBMLMap SBMLid="GAP" COPASIkey="Metabolite_25"/>
    <SBMLMap SBMLid="GDH" COPASIkey="Reaction_4"/>
    <SBMLMap SBMLid="GL6P" COPASIkey="Metabolite_27"/>
    <SBMLMap SBMLid="GL6P_HYDROLYSIS" COPASIkey="Reaction_67"/>
    <SBMLMap SBMLid="GLC_feed" COPASIkey="Reaction_55"/>
    <SBMLMap SBMLid="GLCp" COPASIkey="Metabolite_149"/>
    <SBMLMap SBMLid="GLCx" COPASIkey="Metabolite_141"/>
    <SBMLMap SBMLid="GLT" COPASIkey="Reaction_25"/>
    <SBMLMap SBMLid="GLX" COPASIkey="Metabolite_29"/>
    <SBMLMap SBMLid="GND" COPASIkey="Reaction_11"/>
    <SBMLMap SBMLid="GPM" COPASIkey="Reaction_6"/>
    <SBMLMap SBMLid="GROWTH" COPASIkey="Reaction_59"/>
    <SBMLMap SBMLid="H2O" COPASIkey="Metabolite_131"/>
    <SBMLMap SBMLid="HCO3" COPASIkey="Metabolite_111"/>
    <SBMLMap SBMLid="Hin" COPASIkey="Metabolite_129"/>
    <SBMLMap SBMLid="Hout" COPASIkey="Metabolite_147"/>
    <SBMLMap SBMLid="ICD" COPASIkey="Reaction_28"/>
    <SBMLMap SBMLid="ICIT" COPASIkey="Metabolite_31"/>
    <SBMLMap SBMLid="KDPG" COPASIkey="Metabolite_33"/>
    <SBMLMap SBMLid="KdADPMg" COPASIkey="ModelValue_1"/>
    <SBMLMap SBMLid="KdATPMg" COPASIkey="ModelValue_2"/>
    <SBMLMap SBMLid="KdFDPMg" COPASIkey="ModelValue_3"/>
    <SBMLMap SBMLid="KeqNDH" COPASIkey="ModelValue_7"/>
    <SBMLMap SBMLid="KmACO_ACN" COPASIkey="ModelValue_6"/>
    <SBMLMap SBMLid="KmCIT_ACN" COPASIkey="ModelValue_5"/>
    <SBMLMap SBMLid="KmICIT_ACN" COPASIkey="ModelValue_4"/>
    <SBMLMap SBMLid="LPD" COPASIkey="Reaction_29"/>
    <SBMLMap SBMLid="MAD" COPASIkey="Reaction_23"/>
    <SBMLMap SBMLid="MAL" COPASIkey="Metabolite_35"/>
    <SBMLMap SBMLid="MDH" COPASIkey="Reaction_34"/>
    <SBMLMap SBMLid="MG" COPASIkey="Metabolite_115"/>
    <SBMLMap SBMLid="MN" COPASIkey="Metabolite_127"/>
    <SBMLMap SBMLid="MQO" COPASIkey="Reaction_33"/>
    <SBMLMap SBMLid="MgADP" COPASIkey="Metabolite_117"/>
    <SBMLMap SBMLid="MgATP" COPASIkey="Metabolite_119"/>
    <SBMLMap SBMLid="MgFDP" COPASIkey="Metabolite_121"/>
    <SBMLMap SBMLid="NAD" COPASIkey="Metabolite_37"/>
    <SBMLMap SBMLid="NADH" COPASIkey="Metabolite_39"/>
    <SBMLMap SBMLid="NADH_req" COPASIkey="Reaction_41"/>
    <SBMLMap SBMLid="NADP" COPASIkey="Metabolite_41"/>
    <SBMLMap SBMLid="NADPH" COPASIkey="Metabolite_43"/>
    <SBMLMap SBMLid="NDHII" COPASIkey="Reaction_58"/>
    <SBMLMap SBMLid="O2" COPASIkey="Metabolite_133"/>
    <SBMLMap SBMLid="OAA" COPASIkey="Metabolite_45"/>
    <SBMLMap SBMLid="P" COPASIkey="Metabolite_113"/>
    <SBMLMap SBMLid="PCK" COPASIkey="Reaction_21"/>
    <SBMLMap SBMLid="PDH" COPASIkey="Reaction_24"/>
    <SBMLMap SBMLid="PEP" COPASIkey="Metabolite_47"/>
    <SBMLMap SBMLid="PFK" COPASIkey="Reaction_1"/>
    <SBMLMap SBMLid="PGA2" COPASIkey="Metabolite_49"/>
    <SBMLMap SBMLid="PGA3" COPASIkey="Metabolite_51"/>
    <SBMLMap SBMLid="PGI" COPASIkey="Reaction_0"/>
    <SBMLMap SBMLid="PGK" COPASIkey="Reaction_5"/>
    <SBMLMap SBMLid="PGL" COPASIkey="Reaction_10"/>
    <SBMLMap SBMLid="PGN" COPASIkey="Metabolite_53"/>
    <SBMLMap SBMLid="PIT" COPASIkey="Reaction_62"/>
    <SBMLMap SBMLid="PNT_req" COPASIkey="Reaction_42"/>
    <SBMLMap SBMLid="PPC" COPASIkey="Reaction_20"/>
    <SBMLMap SBMLid="PPS" COPASIkey="Reaction_22"/>
    <SBMLMap SBMLid="PTA" COPASIkey="Reaction_49"/>
    <SBMLMap SBMLid="PTS_0" COPASIkey="Reaction_50"/>
    <SBMLMap SBMLid="PTS_1" COPASIkey="Reaction_51"/>
    <SBMLMap SBMLid="PTS_2" COPASIkey="Reaction_52"/>
    <SBMLMap SBMLid="PTS_3" COPASIkey="Reaction_53"/>
    <SBMLMap SBMLid="PTS_4" COPASIkey="Reaction_54"/>
    <SBMLMap SBMLid="PYK" COPASIkey="Reaction_8"/>
    <SBMLMap SBMLid="PYR" COPASIkey="Metabolite_55"/>
    <SBMLMap SBMLid="Pp" COPASIkey="Metabolite_151"/>
    <SBMLMap SBMLid="Px" COPASIkey="Metabolite_143"/>
    <SBMLMap SBMLid="Q" COPASIkey="Metabolite_57"/>
    <SBMLMap SBMLid="QH2" COPASIkey="Metabolite_59"/>
    <SBMLMap SBMLid="R5P" COPASIkey="Metabolite_61"/>
    <SBMLMap SBMLid="RPE" COPASIkey="Reaction_12"/>
    <SBMLMap SBMLid="RPI" COPASIkey="Reaction_13"/>
    <SBMLMap SBMLid="RU5P" COPASIkey="Metabolite_63"/>
    <SBMLMap SBMLid="S7P" COPASIkey="Metabolite_65"/>
    <SBMLMap SBMLid="S7P_E4P_TAL" COPASIkey="Reaction_18"/>
    <SBMLMap SBMLid="S7P_R5P_TKT" COPASIkey="Reaction_16"/>
    <SBMLMap SBMLid="SDH" COPASIkey="Reaction_31"/>
    <SBMLMap SBMLid="SK" COPASIkey="Reaction_30"/>
    <SBMLMap SBMLid="SQR" COPASIkey="Reaction_57"/>
    <SBMLMap SBMLid="SUC" COPASIkey="Metabolite_67"/>
    <SBMLMap SBMLid="SUCCOA" COPASIkey="Metabolite_69"/>
    <SBMLMap SBMLid="TPI" COPASIkey="Reaction_3"/>
    <SBMLMap SBMLid="X5P" COPASIkey="Metabolite_71"/>
    <SBMLMap SBMLid="X5P_GAP_TKT" COPASIkey="Reaction_14"/>
    <SBMLMap SBMLid="XCH_ACE1" COPASIkey="Reaction_64"/>
    <SBMLMap SBMLid="XCH_ACE2" COPASIkey="Reaction_66"/>
    <SBMLMap SBMLid="XCH_GLC" COPASIkey="Reaction_61"/>
    <SBMLMap SBMLid="XCH_P" COPASIkey="Reaction_63"/>
    <SBMLMap SBMLid="ZWF" COPASIkey="Reaction_9"/>
    <SBMLMap SBMLid="_ACE_OUT" COPASIkey="Reaction_65"/>
    <SBMLMap SBMLid="cell" COPASIkey="Compartment_1"/>
    <SBMLMap SBMLid="cell_periplasm" COPASIkey="Compartment_5"/>
    <SBMLMap SBMLid="ei" COPASIkey="Metabolite_73"/>
    <SBMLMap SBMLid="eiP" COPASIkey="Metabolite_75"/>
    <SBMLMap SBMLid="eiia" COPASIkey="Metabolite_77"/>
    <SBMLMap SBMLid="eiiaP" COPASIkey="Metabolite_79"/>
    <SBMLMap SBMLid="eiicb" COPASIkey="Metabolite_81"/>
    <SBMLMap SBMLid="eiicbP" COPASIkey="Metabolite_83"/>
    <SBMLMap SBMLid="extracellular" COPASIkey="Compartment_3"/>
    <SBMLMap SBMLid="hpr" COPASIkey="Metabolite_85"/>
    <SBMLMap SBMLid="hprP" COPASIkey="Metabolite_87"/>
    <SBMLMap SBMLid="icd" COPASIkey="Metabolite_89"/>
    <SBMLMap SBMLid="icdP" COPASIkey="Metabolite_91"/>
    <SBMLMap SBMLid="tal" COPASIkey="Metabolite_93"/>
    <SBMLMap SBMLid="talC3" COPASIkey="Metabolite_95"/>
    <SBMLMap SBMLid="tkt" COPASIkey="Metabolite_97"/>
    <SBMLMap SBMLid="tktC2" COPASIkey="Metabolite_99"/>
  </SBMLReference>
</COPASI>
