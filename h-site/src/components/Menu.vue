<template>
 <div @mouseover="active=true" @mouseleave="active=false" id="full">
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
				},
        {
          id: 4,
          name: "Profile",
          htmlRef: "/profile"
        },
        {
          id: 5,
          name: "Blogg",
          htmlRef: "/blog"
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
#full {
  width: 1rem;
  height: 1rem;
  position: relative;
}
#contents {
  position: relative;
  left: 2px;
  top: 5px;
  font-size: 1em;
  height: var(--nav-height);
  background-color: white;
  width: 165px;
  border-radius: 8px;
  }
#menu {
  position: relative;
  top: 5px;
  font-size: 1.25em;
  padding: 9px;
  height: var(--nav-height);
  border-radius: 8px;
  background-color: black;
  width: 90px;
  color: blue;
}
ul {
  text-align: left;
  text-decoration: none;
  color: black;
  background-color: snow;
  list-style-type: none;
  font-family: Trebuchet MS, serif;
  }
a {
  color: red;
  text-decoration: none;
  cursor: pointer;
}
</style>