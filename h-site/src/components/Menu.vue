<template>
 <div @mouseover="active=true" @mouseleave="active=false">
  <div id="contents-button" v-if="!active">
   <button id="menu">Menu</button>
  </div>
  <div id="contents" v-if="active">
   <ul id="content-list"> 
    <li v-for="(element, idx) in elements" v-bind:key="idx"> 
      <a @click="changeRoute(element.htmlRef)" class="content-links"> {{ element.name }} </a> 
    </li>
   </ul> 
   </div>
 </div>
</template>

<script>
export default {
	data() {
		return {
			elements: [
				{
					id: 1,
					name: "Home",
					htmlRef: "/"
				},
				{
					id: 2,
					name: "Todo Manager",
					htmlRef: "/todos"
				},
				{
					id: 3,
					name: "Quiz",
					htmlRef: "/quiz"
				}],
			active: false		
		}
	},
	methods: {
		stuff() {
			console.log("Stuff")
			this.active = !this.active
		},
    changeRoute(newRoute) {
      this.$router.push(newRoute)
    }
	}
}
</script>

<style>
#contents {
  position: relative;
  left: 2px;
  top: 5px;
  font-size: 1em;
  height: var(--nav-height);
  background-color: lightgrey;
  width: 165px;
  border-radius: 8px;
  }
#menu {
  position: relative;
  left: -10rem;
  top: 5px;
  font-size: 1.25em;
  padding: 9px;
  font-family: Trebuchet MS, serif;
  height: var(--nav-height);
  border-radius: 8px;
  background-color: lightgrey;
  width: 90px;
}
ul {
  text-align: left;
  text-decoration: none;
  color: black;
  background-color: lightgrey;
  list-style-type: none;
  font-family: Trebuchet MS, serif;
  }
a {
  color: red;
  text-decoration: none;
  cursor: pointer;
}
</style>