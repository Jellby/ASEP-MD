<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns="http://www.w3.org/1999/xhtml">

<xsl:include href="xhtml.xsl"/>

<xsl:variable name="title">
  <xsl:copy-of select="document/title[@lang=$lang]"/>
</xsl:variable>

<xsl:template match="document">
  <xsl:element name="html">
    <xsl:call-template name="head"/>

    <xsl:element name="body">

      <xsl:element name="h2"><xsl:value-of select="$title"/></xsl:element>

      <xsl:apply-templates select="paragraph[@lang=$lang] | sample[@lang=$lang]"/>

      <xsl:apply-templates select="quantity"/>

    </xsl:element>
  </xsl:element>
</xsl:template>

<xsl:template match="sample">
  <xsl:element name="div">
    <xsl:attribute name="class">sample</xsl:attribute>
    <xsl:apply-templates/>
  </xsl:element>
</xsl:template>

<xsl:template match="description | paragraph">
  <xsl:element name="p">
    <xsl:attribute name="class">description</xsl:attribute>
    <xsl:apply-templates/>
  </xsl:element>
</xsl:template>

<xsl:template match="quantity">
  <xsl:element name="div">
    <xsl:attribute name="class">quantity</xsl:attribute>

    <xsl:apply-templates select="name[@lang=$lang]"/>

    <xsl:apply-templates select="description[@lang=$lang] | description[not(@lang)]"/>

  </xsl:element>
</xsl:template>

<xsl:template match="name">
  <xsl:element name="div">
    <xsl:attribute name="class">names</xsl:attribute>
    <xsl:call-template name="var">
      <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
    </xsl:call-template>
  </xsl:element>
</xsl:template>

</xsl:stylesheet>
