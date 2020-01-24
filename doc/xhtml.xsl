<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns="http://www.w3.org/1999/xhtml">

<xsl:include href="functions.xsl"/>

<xsl:output method="xml" indent="yes" doctype-system="http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd" doctype-public="-//W3C//DTD XHTML 1.1//EN" encoding="utf-8"/>

<xsl:template name="head">
  <xsl:element name="head">
    <xsl:element name="title"><xsl:value-of select="$title"/></xsl:element>
    <xsl:element name="link">
      <xsl:attribute name="rel">stylesheet</xsl:attribute>
      <xsl:attribute name="type">text/css</xsl:attribute>
      <xsl:attribute name="href">style.css</xsl:attribute>
    </xsl:element>
    <xsl:element name="script">
      <xsl:attribute name="type">text/x-mathjax-config</xsl:attribute>
MathJax.Hub.Config({
  "HTML-CSS": {
    mtextFontInherit: true
  }
});
    </xsl:element>
    <xsl:element name="script">
      <xsl:attribute name="type">text/javascript</xsl:attribute>
      <xsl:attribute name="src">http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML</xsl:attribute>
    </xsl:element>
  </xsl:element>
</xsl:template>

<xsl:template match="math">
  <xsl:text>\(</xsl:text>
  <xsl:choose>
    <xsl:when test="@latex">
      <xsl:value-of select="@latex"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="."/>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>\)</xsl:text>
</xsl:template>

<xsl:template match="math" mode="equation">
  <xsl:text>\(\displaystyle </xsl:text>
  <xsl:choose>
    <xsl:when test="@latex">
      <xsl:value-of select="@latex"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="."/>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>\)</xsl:text>
</xsl:template>

<xsl:template match="exp">
  <xsl:text>·10</xsl:text>
  <xsl:element name="sup"><xsl:apply-templates/></xsl:element>
</xsl:template>

<xsl:template name="tag">
  <xsl:param name="text"/>
  <xsl:element name="i">
    <xsl:attribute name="class">tag</xsl:attribute>
    <xsl:copy-of select="$text"/>
  </xsl:element>
</xsl:template>

<xsl:template name="val">
  <xsl:param name="text"/>
  <xsl:element name="tt">
    <xsl:attribute name="class">val</xsl:attribute>
    <xsl:copy-of select="$text"/>
  </xsl:element>
</xsl:template>
<xsl:template match="val">
  <xsl:call-template name="val">
    <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
  </xsl:call-template>
</xsl:template>

<xsl:template name="var">
  <xsl:param name="text"/>
  <xsl:element name="b">
    <xsl:attribute name="class">var</xsl:attribute>
    <xsl:copy-of select="$text"/>
  </xsl:element>
</xsl:template>
<xsl:template match="var">
  <xsl:call-template name="var">
    <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
  </xsl:call-template>
</xsl:template>

<xsl:template name="prog">
  <xsl:param name="text"/>
  <xsl:element name="b">
    <xsl:attribute name="class">prog</xsl:attribute>
    <xsl:copy-of select="$text"/>
  </xsl:element>
</xsl:template>
<xsl:template match="prog">
  <xsl:call-template name="prog">
    <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
  </xsl:call-template>
</xsl:template>

<xsl:template match="equation">
  <xsl:element name="span">
    <xsl:attribute name="class">equation</xsl:attribute>
    <xsl:if test="@id">
      <xsl:attribute name="id"><xsl:value-of select="@id"/></xsl:attribute>
      <xsl:attribute name="title"><xsl:value-of select="@id"/></xsl:attribute>
      <xsl:element name="span">
        <xsl:attribute name="class">label</xsl:attribute>
        <xsl:text>(</xsl:text>
        <xsl:call-template name="eqlabel">
          <xsl:with-param name="label" select="@id"/>
        </xsl:call-template>
        <xsl:text>)</xsl:text>
      </xsl:element>
    </xsl:if>
    <xsl:apply-templates mode="equation"/>
  </xsl:element>
</xsl:template>

<xsl:template match="ref">
  <xsl:element name="a">
    <xsl:attribute name="class">eqref</xsl:attribute>
    <xsl:attribute name="href"><xsl:text>#</xsl:text><xsl:value-of select="@label"/></xsl:attribute>
    <xsl:text>(</xsl:text>
    <xsl:call-template name="eqlabel">
      <xsl:with-param name="label" select="@label"/>
    </xsl:call-template>
    <xsl:text>)</xsl:text>
  </xsl:element>
</xsl:template>

</xsl:stylesheet>
