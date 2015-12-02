/**
 * Copyright (C) 2015 - present by OpenGamma Inc. and the OpenGamma group of companies
 * <p>
 * Please see distribution for license.
 */
package com.opengamma.strata.market.key;

import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.joda.beans.Bean;
import org.joda.beans.BeanDefinition;
import org.joda.beans.ImmutableBean;
import org.joda.beans.JodaBeanUtils;
import org.joda.beans.MetaProperty;
import org.joda.beans.Property;
import org.joda.beans.PropertyDefinition;
import org.joda.beans.impl.direct.DirectFieldsBeanBuilder;
import org.joda.beans.impl.direct.DirectMetaBean;
import org.joda.beans.impl.direct.DirectMetaProperty;
import org.joda.beans.impl.direct.DirectMetaPropertyMap;

import com.opengamma.strata.basics.market.FieldName;
import com.opengamma.strata.basics.market.MarketDataBox;
import com.opengamma.strata.basics.market.ScenarioMarketDataKey;
import com.opengamma.strata.collect.array.DoubleArray;
import com.opengamma.strata.collect.id.StandardId;
import com.opengamma.strata.market.value.QuotesArray;

/**
 * A key identifying a {@link QuotesArray} containing values for a piece of quoted market data in multiple scenarios.
 */
@BeanDefinition
public final class QuotesArrayKey
    implements ScenarioMarketDataKey<Double, QuotesArray>, ImmutableBean {

  /**
   * The ID of the market data that is required, typically an ID from an external data provider.
   * For example, 'Bloomberg~AAPL'.
   */
  @PropertyDefinition(validate = "notNull")
  private final StandardId id;
  /**
   * The field name in the market data record that is required.
   * For example, {@link FieldName#MARKET_VALUE}.
   */
  @PropertyDefinition(validate = "notNull")
  private final FieldName fieldName;

  /**
   * Returns a key identifying the market data with the specified ID and field name.
   *
   * @return a key identifying the market data with the specified ID and field name
   */
  public static QuotesArrayKey of(StandardId id, FieldName fieldName) {
    return new QuotesArrayKey(id, fieldName);
  }

  @Override
  public QuoteKey getMarketDataKey() {
    return QuoteKey.of(id, fieldName);
  }

  @Override
  public Class<QuotesArray> getScenarioMarketDataType() {
    return QuotesArray.class;
  }

  @Override
  public QuotesArray createScenarioValue(MarketDataBox<Double> marketDataBox) {
    // TODO This suggests MarketDataBox needs a stream() method or at least a helper to create a stream from a box
    double[] quotes = new double[marketDataBox.getScenarioCount()];

    for (int i = 0; i < marketDataBox.getScenarioCount(); i++) {
      quotes[i] = marketDataBox.getValue(i);
    }
    return QuotesArray.of(DoubleArray.ofUnsafe(quotes));
  }

  //------------------------- AUTOGENERATED START -------------------------
  ///CLOVER:OFF
  /**
   * The meta-bean for {@code QuotesArrayKey}.
   * @return the meta-bean, not null
   */
  public static QuotesArrayKey.Meta meta() {
    return QuotesArrayKey.Meta.INSTANCE;
  }

  static {
    JodaBeanUtils.registerMetaBean(QuotesArrayKey.Meta.INSTANCE);
  }

  /**
   * Returns a builder used to create an instance of the bean.
   * @return the builder, not null
   */
  public static QuotesArrayKey.Builder builder() {
    return new QuotesArrayKey.Builder();
  }

  private QuotesArrayKey(
      StandardId id,
      FieldName fieldName) {
    JodaBeanUtils.notNull(id, "id");
    JodaBeanUtils.notNull(fieldName, "fieldName");
    this.id = id;
    this.fieldName = fieldName;
  }

  @Override
  public QuotesArrayKey.Meta metaBean() {
    return QuotesArrayKey.Meta.INSTANCE;
  }

  @Override
  public <R> Property<R> property(String propertyName) {
    return metaBean().<R>metaProperty(propertyName).createProperty(this);
  }

  @Override
  public Set<String> propertyNames() {
    return metaBean().metaPropertyMap().keySet();
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the ID of the market data that is required, typically an ID from an external data provider.
   * For example, 'Bloomberg~AAPL'.
   * @return the value of the property, not null
   */
  public StandardId getId() {
    return id;
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the field name in the market data record that is required.
   * For example, {@link FieldName#MARKET_VALUE}.
   * @return the value of the property, not null
   */
  public FieldName getFieldName() {
    return fieldName;
  }

  //-----------------------------------------------------------------------
  /**
   * Returns a builder that allows this bean to be mutated.
   * @return the mutable builder, not null
   */
  public Builder toBuilder() {
    return new Builder(this);
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == this) {
      return true;
    }
    if (obj != null && obj.getClass() == this.getClass()) {
      QuotesArrayKey other = (QuotesArrayKey) obj;
      return JodaBeanUtils.equal(id, other.id) &&
          JodaBeanUtils.equal(fieldName, other.fieldName);
    }
    return false;
  }

  @Override
  public int hashCode() {
    int hash = getClass().hashCode();
    hash = hash * 31 + JodaBeanUtils.hashCode(id);
    hash = hash * 31 + JodaBeanUtils.hashCode(fieldName);
    return hash;
  }

  @Override
  public String toString() {
    StringBuilder buf = new StringBuilder(96);
    buf.append("QuotesArrayKey{");
    buf.append("id").append('=').append(id).append(',').append(' ');
    buf.append("fieldName").append('=').append(JodaBeanUtils.toString(fieldName));
    buf.append('}');
    return buf.toString();
  }

  //-----------------------------------------------------------------------
  /**
   * The meta-bean for {@code QuotesArrayKey}.
   */
  public static final class Meta extends DirectMetaBean {
    /**
     * The singleton instance of the meta-bean.
     */
    static final Meta INSTANCE = new Meta();

    /**
     * The meta-property for the {@code id} property.
     */
    private final MetaProperty<StandardId> id = DirectMetaProperty.ofImmutable(
        this, "id", QuotesArrayKey.class, StandardId.class);
    /**
     * The meta-property for the {@code fieldName} property.
     */
    private final MetaProperty<FieldName> fieldName = DirectMetaProperty.ofImmutable(
        this, "fieldName", QuotesArrayKey.class, FieldName.class);
    /**
     * The meta-properties.
     */
    private final Map<String, MetaProperty<?>> metaPropertyMap$ = new DirectMetaPropertyMap(
        this, null,
        "id",
        "fieldName");

    /**
     * Restricted constructor.
     */
    private Meta() {
    }

    @Override
    protected MetaProperty<?> metaPropertyGet(String propertyName) {
      switch (propertyName.hashCode()) {
        case 3355:  // id
          return id;
        case 1265009317:  // fieldName
          return fieldName;
      }
      return super.metaPropertyGet(propertyName);
    }

    @Override
    public QuotesArrayKey.Builder builder() {
      return new QuotesArrayKey.Builder();
    }

    @Override
    public Class<? extends QuotesArrayKey> beanType() {
      return QuotesArrayKey.class;
    }

    @Override
    public Map<String, MetaProperty<?>> metaPropertyMap() {
      return metaPropertyMap$;
    }

    //-----------------------------------------------------------------------
    /**
     * The meta-property for the {@code id} property.
     * @return the meta-property, not null
     */
    public MetaProperty<StandardId> id() {
      return id;
    }

    /**
     * The meta-property for the {@code fieldName} property.
     * @return the meta-property, not null
     */
    public MetaProperty<FieldName> fieldName() {
      return fieldName;
    }

    //-----------------------------------------------------------------------
    @Override
    protected Object propertyGet(Bean bean, String propertyName, boolean quiet) {
      switch (propertyName.hashCode()) {
        case 3355:  // id
          return ((QuotesArrayKey) bean).getId();
        case 1265009317:  // fieldName
          return ((QuotesArrayKey) bean).getFieldName();
      }
      return super.propertyGet(bean, propertyName, quiet);
    }

    @Override
    protected void propertySet(Bean bean, String propertyName, Object newValue, boolean quiet) {
      metaProperty(propertyName);
      if (quiet) {
        return;
      }
      throw new UnsupportedOperationException("Property cannot be written: " + propertyName);
    }

  }

  //-----------------------------------------------------------------------
  /**
   * The bean-builder for {@code QuotesArrayKey}.
   */
  public static final class Builder extends DirectFieldsBeanBuilder<QuotesArrayKey> {

    private StandardId id;
    private FieldName fieldName;

    /**
     * Restricted constructor.
     */
    private Builder() {
    }

    /**
     * Restricted copy constructor.
     * @param beanToCopy  the bean to copy from, not null
     */
    private Builder(QuotesArrayKey beanToCopy) {
      this.id = beanToCopy.getId();
      this.fieldName = beanToCopy.getFieldName();
    }

    //-----------------------------------------------------------------------
    @Override
    public Object get(String propertyName) {
      switch (propertyName.hashCode()) {
        case 3355:  // id
          return id;
        case 1265009317:  // fieldName
          return fieldName;
        default:
          throw new NoSuchElementException("Unknown property: " + propertyName);
      }
    }

    @Override
    public Builder set(String propertyName, Object newValue) {
      switch (propertyName.hashCode()) {
        case 3355:  // id
          this.id = (StandardId) newValue;
          break;
        case 1265009317:  // fieldName
          this.fieldName = (FieldName) newValue;
          break;
        default:
          throw new NoSuchElementException("Unknown property: " + propertyName);
      }
      return this;
    }

    @Override
    public Builder set(MetaProperty<?> property, Object value) {
      super.set(property, value);
      return this;
    }

    @Override
    public Builder setString(String propertyName, String value) {
      setString(meta().metaProperty(propertyName), value);
      return this;
    }

    @Override
    public Builder setString(MetaProperty<?> property, String value) {
      super.setString(property, value);
      return this;
    }

    @Override
    public Builder setAll(Map<String, ? extends Object> propertyValueMap) {
      super.setAll(propertyValueMap);
      return this;
    }

    @Override
    public QuotesArrayKey build() {
      return new QuotesArrayKey(
          id,
          fieldName);
    }

    //-----------------------------------------------------------------------
    /**
     * Sets the ID of the market data that is required, typically an ID from an external data provider.
     * For example, 'Bloomberg~AAPL'.
     * @param id  the new value, not null
     * @return this, for chaining, not null
     */
    public Builder id(StandardId id) {
      JodaBeanUtils.notNull(id, "id");
      this.id = id;
      return this;
    }

    /**
     * Sets the field name in the market data record that is required.
     * For example, {@link FieldName#MARKET_VALUE}.
     * @param fieldName  the new value, not null
     * @return this, for chaining, not null
     */
    public Builder fieldName(FieldName fieldName) {
      JodaBeanUtils.notNull(fieldName, "fieldName");
      this.fieldName = fieldName;
      return this;
    }

    //-----------------------------------------------------------------------
    @Override
    public String toString() {
      StringBuilder buf = new StringBuilder(96);
      buf.append("QuotesArrayKey.Builder{");
      buf.append("id").append('=').append(JodaBeanUtils.toString(id)).append(',').append(' ');
      buf.append("fieldName").append('=').append(JodaBeanUtils.toString(fieldName));
      buf.append('}');
      return buf.toString();
    }

  }

  ///CLOVER:ON
  //-------------------------- AUTOGENERATED END --------------------------
}
