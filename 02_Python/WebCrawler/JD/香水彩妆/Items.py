# -*- coding: utf-8 -*-

# Define here the models for your scraped items
#
# See documentation in:
# http://doc.scrapy.org/en/latest/topics/items.html

from scrapy import Item, Field


class CategoriesItem(Item):
	category = Field()  #分类名称
	url = Field()  #分类url
	class_id= Field() 
	_class=Field() #分类的大类


class ProductsItem(Item):
	name = Field()  #产品名称
	items_url = Field()  #产品url
	class_url=Field()
	_id = Field()  #产品id
	category = Field()  #产品分类
	_class=Field()
	class_id=Field()
	items_id=Field()
	reallyPrice = Field()  #产品价格
	originalPrice = Field()  #原价
	description = Field()  #产品描述
	shopId = Field()  #shop id
	favourableDesc1 = Field()  #优惠描述1
	favourableDesc2 = Field()  #优惠描述2
	goodRateShow= Field()
	poorRateShow= Field()
	generalRateShow= Field()
	goodCount= Field()
	generalCount= Field()
	poorCount= Field()
	commentCount= Field()
	defaultGoodCount= Field()
	averageScore= Field()	


class ShopItem(Item):
	_id = Field()  #店铺名称
	name = Field()  #店铺名称
	url1 = Field()  #店铺url1
	url2 = Field()  #店铺url2
	shopId = Field()  #shop id
	venderId = Field()  #vender id

class HotCommentTagItem(Item):
	name= Field() 
	count= Field() 
	_class= Field() 
	category= Field() 
	productId= Field() 


