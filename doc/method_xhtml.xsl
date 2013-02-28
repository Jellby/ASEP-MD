<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:include href="xhtml.xsl"/>

<xsl:variable name="title">
  <xsl:copy-of select="document/title[@lang=$lang]"/>
</xsl:variable>

<xsl:template match="document">
  <xsl:element name="html">
    <xsl:attribute name="xmlns">http://www.w3.org/1999/xhtml</xsl:attribute>
    <xsl:call-template name="head"/>

    <xsl:element name="body">

      <xsl:element name="h2"><xsl:value-of select="$title"/></xsl:element>

      <xsl:apply-templates select="paragraph[@lang=$lang]"/>
    </xsl:element>
  </xsl:element>
</xsl:template>

<xsl:template match="paragraph">
  <xsl:element name="p">
    <xsl:attribute name="class">description</xsl:attribute>
    <xsl:apply-templates/>
  </xsl:element>
</xsl:template>

</xsl:stylesheet>
