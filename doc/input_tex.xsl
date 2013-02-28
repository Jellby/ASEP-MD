<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:include href="tex.xsl"/>

<xsl:variable name="title">
  <xsl:copy-of select="document/title[@lang=$lang]"/>
</xsl:variable>

<xsl:template match="document">
  <xsl:copy-of select="$begindocument"/>

  <xsl:text>\chapter*{</xsl:text>
  <xsl:copy-of select="$title"/>
  <xsl:text>}</xsl:text>

  <xsl:apply-templates select="description[@lang=$lang]"/>
  <xsl:apply-templates select="variable[name[@lang='es']='Inicio']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='MaxIter']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='InicioVacio']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='ProgramaMM']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='EjecutableMM']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='ProgramaQM']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='EjecutableQM']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='EntradaMM']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='SalidaMM']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='TrayectoriaMM']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='EntradaQM']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='SalidaQM']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='SalidaOpt']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='FChkGaussian']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='DumpextMoldy']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='CargasExternas']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='TipoCargas']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='FicheroCargas']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='NumConfig']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='DistCargas']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='TipoReduccion']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='CorteRF']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='Dielectrica']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='TipoCavidad']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='RadioCavidad']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='Subdivisiones']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='RadioDisolvente']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='TipoMalla']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='FactorMalla']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='DistMalla']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='OptContinuacion']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='TipoCoordenadas']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='AltCoordenadas']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='MaxIterOpt']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='MetodoOptim']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='HessInicial']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='CalcHessiana']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='Actualizacion']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='BusquedaLineal']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='EstadoTransicion']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='ConvGradOpt']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='ConvPasoOpt']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='ConvEnerOpt']"/>
  <xsl:apply-templates select="variable[name[@lang='es']='MaxPasoOpt']"/>

  <xsl:copy-of select="$enddocument"/>
</xsl:template>

<xsl:template match="variable">
  <xsl:value-of select="$line"/>
  <xsl:value-of select="$newline"/>
  <xsl:text>\noindent\begin{minipage}{\columnwidth}</xsl:text>
  <xsl:value-of select="$newline"/>
  <xsl:call-template name="names"/>

  <xsl:apply-templates select="default[@lang=$lang] | default[not(@lang)]"/>

  <xsl:apply-templates select="description[@lang=$lang]"/>

  <xsl:if test="value">
    <xsl:call-template name="values"/>
  </xsl:if>
  <xsl:text>\end{minipage}</xsl:text>
  <xsl:value-of select="$newline"/>
  <xsl:value-of select="$newline"/>
</xsl:template>

<xsl:template name="names">
  <xsl:apply-templates select="name">
    <xsl:sort select="@lang"/>
  </xsl:apply-templates>
</xsl:template>

<xsl:template match="name">
  <xsl:text>\noindent </xsl:text>
  <xsl:call-template name="tag">
    <xsl:with-param name="text"><xsl:value-of select="@lang"/></xsl:with-param>
  </xsl:call-template>
  <xsl:text>: </xsl:text>
  <xsl:call-template name="var">
    <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
  </xsl:call-template>
  <xsl:text>\\</xsl:text>
  <xsl:value-of select="$newline"/>
</xsl:template>

<xsl:template match="default">
  <xsl:text>\hspace*{3em}(</xsl:text>
  <xsl:call-template name="tag">
    <xsl:with-param name="text"><xsl:value-of select="$tag_default"/></xsl:with-param>
  </xsl:call-template>
  <xsl:text>: </xsl:text>
  <xsl:choose>
    <xsl:when test="not(string(.))">
      <xsl:text>\hspace{3em}</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:choose>
        <xsl:when test="child::exp">
          <xsl:call-template name="val">
            <xsl:with-param name="text"><xsl:text>$</xsl:text><xsl:apply-templates/><xsl:text>$</xsl:text></xsl:with-param>
          </xsl:call-template>
        </xsl:when>
        <xsl:otherwise>
          <xsl:call-template name="val">
            <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
          </xsl:call-template>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>)</xsl:text>
  <xsl:value-of select="$newline"/>
  <xsl:value-of select="$newline"/>
</xsl:template>

<xsl:template match="description">
  <xsl:apply-templates/>
  <xsl:value-of select="$newline"/>
  <xsl:value-of select="$newline"/>
</xsl:template>

<xsl:template match="description" mode="list">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template name="values">
  <xsl:text>\medskip</xsl:text>
  <xsl:value-of select="$newline"/>
  <xsl:call-template name="tag">
    <xsl:with-param name="text"><xsl:value-of select="$tag_allowed"/></xsl:with-param>
  </xsl:call-template>
  <xsl:text>:</xsl:text>
  <xsl:value-of select="$newline"/>
  <xsl:text>\begin{itemize}</xsl:text>
  <xsl:value-of select="$newline"/>
  <xsl:apply-templates select="value"/>
  <xsl:text>\end{itemize}</xsl:text>
  <xsl:value-of select="$newline"/>
  <xsl:value-of select="$newline"/>
</xsl:template>

<xsl:template match="value">
  <xsl:text>\item </xsl:text>
  <xsl:apply-templates select="option"/>
  <xsl:text>: </xsl:text>
  <xsl:apply-templates select="description[@lang=$lang]" mode="list"/>
  <xsl:value-of select="$newline"/>
</xsl:template>

<xsl:template match="option">
  <xsl:choose>
    <xsl:when test="not(string(.))">
      <xsl:call-template name="tag">
        <xsl:with-param name="text"><xsl:value-of select="$tag_other"/></xsl:with-param>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="val">
        <xsl:with-param name="text"><xsl:apply-templates/></xsl:with-param>
      </xsl:call-template>
      <xsl:if test="position() != last()"><xsl:text>, </xsl:text></xsl:if>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

</xsl:stylesheet>
