/**
 * Copyright (C) 2015 - present by OpenGamma Inc. and the OpenGamma group of companies
 * <p>
 * Please see distribution for license.
 */
package com.opengamma.strata.market.key;

import static org.assertj.core.api.Assertions.assertThat;

import org.testng.annotations.Test;

import com.opengamma.strata.basics.market.FieldName;
import com.opengamma.strata.basics.market.MarketDataBox;
import com.opengamma.strata.collect.array.DoubleArray;
import com.opengamma.strata.collect.id.StandardId;
import com.opengamma.strata.market.value.QuotesArray;

@Test
public class QuotesArrayKeyTest {

  private static final QuotesArrayKey KEY = QuotesArrayKey.of(StandardId.of("test", "1"), FieldName.of("fieldName"));

  public void getMarketDataKey() {
    assertThat(KEY.getMarketDataKey()).isEqualTo(QuoteKey.of(StandardId.of("test", "1"), FieldName.of("fieldName")));
  }

  public void getMarketDataType() {
    assertThat(KEY.getScenarioMarketDataType()).isEqualTo(QuotesArray.class);
  }

  public void createScenarioValue() {
    MarketDataBox<Double> box = MarketDataBox.ofScenarioValues(1d, 2d, 3d);
    QuotesArray quotesArray = KEY.createScenarioValue(box);
    assertThat(quotesArray.getValues()).isEqualTo(DoubleArray.of(1d, 2d, 3d));
  }
}
