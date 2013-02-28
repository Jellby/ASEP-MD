<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:include href="functions.xsl"/>

<xsl:output method="text" indent="no"/>

<xsl:variable name="newline">
  <xsl:text>
</xsl:text>
</xsl:variable>

<xsl:variable name="line">
  <xsl:text>================================================================================
</xsl:text>
</xsl:variable>

<xsl:template match="math">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="exp">
  <xsl:text>·10^</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template name="tag">
  <xsl:param name="text"/>
  <xsl:copy-of select="$text"/>
</xsl:template>

<xsl:template name="val">
  <xsl:param name="text"/>
  <xsl:text>¦</xsl:text><xsl:copy-of select="$text"/><xsl:text>¦</xsl:text>
</xsl:template>
<xsl:template match="val">
  <xsl:call-template name="val">
    <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
  </xsl:call-template>
</xsl:template>

<xsl:template name="var">
  <xsl:param name="text"/>
  <xsl:text>¤</xsl:text><xsl:copy-of select="$text"/><xsl:text>¤</xsl:text>
</xsl:template>
<xsl:template match="var">
  <xsl:call-template name="var">
    <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
  </xsl:call-template>
</xsl:template>

<xsl:template name="prog">
  <xsl:param name="text"/>
  <xsl:copy-of select="$text"/>
</xsl:template>
<xsl:template match="prog">
  <xsl:call-template name="prog">
    <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
  </xsl:call-template>
</xsl:template>

<xsl:template match="equation">
  <xsl:text>¬{</xsl:text>
  <xsl:apply-templates/>
  <xsl:if test="@id">
    <xsl:text>¬(</xsl:text>
    <xsl:call-template name="eqlabel">
      <xsl:with-param name="label" select="@id"/>
    </xsl:call-template>
    <xsl:text>)</xsl:text>
  </xsl:if>
  <xsl:text>¬}</xsl:text>
</xsl:template>

<xsl:template match="ref">
  <xsl:text>(</xsl:text>
  <xsl:call-template name="eqlabel">
    <xsl:with-param name="label" select="@label"/>
  </xsl:call-template>
  <xsl:text>)</xsl:text>
</xsl:template>

</xsl:stylesheet>
